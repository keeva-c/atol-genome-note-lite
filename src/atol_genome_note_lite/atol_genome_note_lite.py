#!/usr/bin/env python3

import argparse
from importlib.resources import files
import json
from pathlib import Path
from re import sub
import yaml

import jinja2
import pandas as pd

from atol_genome_note_lite.common import _existing_file, logger
from atol_genome_note_lite.models import (
    BuscoStats,
    SeqStats,
    TaxonomyInfo,
    MerquryStats,
    SoftwareVersions,
)
from canopy_client import CanopySession
from yaml_manifest import Manifest


class IgnoreUndefined(jinja2.Undefined):
    def _fail_with_undefined_error(self, *args, **kwargs):
        logger.warning(self._undefined_message)
        return None


def convert_to_snake_case(str: str) -> str:
    # Defined here:
    # https://linkml.io/linkml/_modules/linkml/utils/helpers.html#convert_to_snake_case,
    # However, the Pydantic generator doesn't convert case so we modify it.
    str = sub(r"(?<=[a-z])(?=[A-Z])|[^a-zA-Z]", " ", str).strip().replace(" ", "_")
    return "".join(str)


def parse_arguments() -> argparse.Namespace:
    my_files = files()

    parser = argparse.ArgumentParser(
        description="Render the templates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # inputs
    _ = parser.add_argument(
        "--template",
        type=_existing_file,
        default=Path(my_files, "templates", "base.md.j2"),
        help="Base template.",
    )

    _ = parser.add_argument(
        "--seq_stats",
        type=_existing_file,
        help="Output from seq stats.",
    )

    _ = parser.add_argument(
        "--busco_hap1",
        type=_existing_file,
        help="short_summary JSON file from BUSCO (hap1)",
    )

    _ = parser.add_argument(
        "--busco_hap2",
        type=_existing_file,
        help="short_summary JSON file from BUSCO (hap2)",
    )

    _ = parser.add_argument(
        "--merqury_completeness",
        type=_existing_file,
        help="asm.ccs.completeness.stats file from Merqury.",
    )

    _ = parser.add_argument(
        "--merqury_qv",
        type=_existing_file,
        help="asm.ccs.qv file from Merqury",
    )

    _ = parser.add_argument(
        "--ga_ver",
        type=_existing_file,
        help="genomeassembly software versions yaml file."
    )

    _ = parser.add_argument(
        "--ascc_ver",
        type=_existing_file,
        help="ascc software versions yaml file."
    )

    _ = parser.add_argument(
        "--treeval_ver",
        type=_existing_file,
        help="treeval software versions yaml file."
    )

    # outputs
    _ = parser.add_argument(
        "--outdir", type=Path, help="Output directory.", required=True
    )

    _ = parser.add_argument(
        "manifest", type=_existing_file, help="AToL manifest for this assembly."
    )

    return parser.parse_args()


def recursive_rekey(d, rename_func=convert_to_snake_case):
    if not isinstance(d, dict):
        return d
    return {
        rename_func(k): (recursive_rekey(v, rename_func) if isinstance(v, dict) else v)
        for k, v in d.items()
    }


def read_merqury_stats(merqury_stats_file: Path) -> dict[str, MerquryStats]:
    merqury_stats = pd.read_csv(merqury_stats_file, sep="\t", index_col="Assembly")
    merqury_stats_dict = merqury_stats.to_dict(orient="index", index=True)
    return {
        k: MerquryStats(**recursive_rekey(v)) for k, v in merqury_stats_dict.items()
    }


def read_seq_stats(seq_stats_file: Path) -> dict[str, SeqStats]:
    """
    Parse the seqkit stats output into a dict of SeqStats models keyed by the
    file path
    """
    seq_stats = pd.read_csv(seq_stats_file, sep="\t", index_col="file")
    seq_stats_dict = seq_stats.to_dict(orient="index", index=True)
    return {k: SeqStats(**v) for k, v in seq_stats_dict.items()}


def read_busco_json(busco_json_file: Path) -> dict[str, BuscoStats]:
    with open(busco_json_file, "rb") as f:
        busco_dict = json.load(f)

    return BuscoStats(**recursive_rekey(busco_dict))


def read_versions_yaml(software_summary: Path, rename_func=convert_to_snake_case) -> SoftwareVersions | None:
    """Parse the software version yaml to extract the pipeline version in the shape of the SoftwareVersions model"""
    with open(software_summary, "rt") as f:
        software = yaml.load(f, Loader=yaml.SafeLoader)
        workflows = software.get("Workflow", {})
        normalised_workflows = {}
        for workflow, version in workflows.items():
            normalised_workflows[rename_func(workflow)] = version

        for pydantic_slot in SoftwareVersions.model_fields:
            if pydantic_slot in normalised_workflows:
                return SoftwareVersions(**{pydantic_slot: normalised_workflows[pydantic_slot]})

        logger.warning(f"Could not find a software version in {software_summary}")
        return None


def render_template_string(template_string: str, **kwargs) -> str:

    env = jinja2.Environment(undefined=IgnoreUndefined)
    template = env.from_string(template_string)
    context = {**kwargs}
    return template.render(context)


def render_template(template_path: Path, **kwargs) -> str:
    loader_path = template_path.parent.as_posix()
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(loader_path), undefined=IgnoreUndefined
    )
    template = env.get_template(template_path.name)
    return template.render(kwargs)


def main():

    args = parse_arguments()
    logger.name = "atol-genome-note-lite"

    # 0. Read the assembly Manifest
    with open(args.manifest, "rb") as f:
        manifest = Manifest.model_validate_json(f.read())
    taxon_id = manifest.taxon_id

    # The context for Assembly stats is a dict with keys hap1 and hap2. The
    # inner dicts have objects for the different result classes, e.g.
    # {"hap1": {
    #     "seq_stats": SeqStats()
    #     "busco_stats": BuscoStats()
    #     ...
    #     },
    # }

    assemblies = {
        "hap1": manifest.treeval_assembly.outputs_for("ascc")
        .get("PRIMARY", Path())
        .as_posix(),
        "hap2": manifest.treeval_assembly.outputs_for("ascc")
        .get("HAPLO", Path())
        .as_posix(),
        "combined": manifest.treeval_assembly.outputs_for("ascc")
        .get("COMBINED", Path())
        .as_posix(),
    }
    assembly_context = {k: {} for k in assemblies.keys()}

    # 1. Read the required information off Canopy
    canopy_session = CanopySession()

    # 1.1. Read taxonomy info off Canopy
    canopy_taxonomy_info = canopy_session.get_taxonomy_info(taxon_id=taxon_id).json()
    taxonomy_info = TaxonomyInfo(**canopy_taxonomy_info)

    # 1.2. Read samples info off Canopy
        # For package in manifest, retrieve sample submission:
            # Take biosample_accession
            # Take sample ID
        # For sample ID in list of sample IDs:
            # Retrieve and store metadata in SampleMetadata class
            # Store biosample_accession in SampleMetadata
            # (not currently relevant) Take derived_from_sample_id values
        # Not currently possible to retrieve the SAMEA accession for specimen-level samples

    # 1.3. Read runs info off Canopy
        # Use assembly ID to retrieve QC reads
            # Store metadata in RunsMetadata class
            # Take experiment IDs

    # 1.4. Read experiment info off Canopy
        # For experiment ID in list of experiment IDs:
            # Retrieve and store metadata in ExperimentMetadata class

    # 1.5. Read assembly info off Canopy
        # Use assembly ID to retrieve assembly metadaata
            # Store metadata in AssemblyMetadata class
    
    # 1.6. Read project info off Canopy
        # Use manifest.taxon_id to retrieve root-level project metadata
            # Store metadata in ProjectMetadata class

    # 2. Harvest info from results files
    if args.seq_stats:
        seq_stats = read_seq_stats(seq_stats_file=args.seq_stats)
        for key, assembly_path in assemblies.items():
            assembly_context[key]["seq_stats"] = seq_stats.get(assembly_path)

    if args.busco_hap1:
        assembly_context["hap1"]["busco_stats"] = read_busco_json(args.busco_hap1)

    if args.busco_hap2:
        assembly_context["hap2"]["busco_stats"] = read_busco_json(args.busco_hap2)

    if args.merqury_completeness:
        merqury_completeness = read_merqury_stats(args.merqury_completeness)
        _hap1_key = [
            key for key in merqury_completeness.keys() if "hap1" in key.lower()
        ][0]
        _hap2_key = [
            key for key in merqury_completeness.keys() if "hap1" in key.lower()
        ][0]
        _combined_key = [
            key
            for key in merqury_completeness.keys()
            if ("+" in key.lower() or "both" in key.lower())
        ][0]

        # jank TODO make this general
        assembly_context["hap1"]["merqury_completeness"] = merqury_completeness[
            _hap1_key
        ]

        assembly_context["hap2"]["merqury_completeness"] = merqury_completeness[
            _hap2_key
        ]

        assembly_context["combined"]["merqury_completeness"] = merqury_completeness[
            _combined_key
        ]

    if args.merqury_qv:
        merqury_qv = read_merqury_stats(args.merqury_qv)
        _hap1_key = [key for key in merqury_qv.keys() if "hap1" in key.lower()][0]
        _hap2_key = [key for key in merqury_qv.keys() if "hap1" in key.lower()][0]
        _combined_key = [
            key
            for key in merqury_qv.keys()
            if ("+" in key.lower() or "both" in key.lower())
        ][0]

        # jank TODO make this general
        assembly_context["hap1"]["merqury_qv"] = merqury_qv[_hap1_key]
        assembly_context["hap2"]["merqury_qv"] = merqury_qv[_hap2_key]
        assembly_context["combined"]["merqury_qv"] = merqury_qv[_combined_key]

    # 2.5. Harvest info from software version files
    method_versions = {}

    if args.ga_ver:
        ga_version = read_versions_yaml(args.ga_ver)
        if ga_version is not None:
            method_versions.update(ga_version.model_dump(exclude_unset=True, exclude_none=True))
    if args.ascc_ver:
        ascc_version = read_versions_yaml(args.ascc_ver)
        if ascc_version is not None:
            method_versions.update(ascc_version.model_dump(exclude_unset=True, exclude_none=True))
    if args.treeval_ver:
        treeval_version = read_versions_yaml(args.treeval_ver)
        if treeval_version is not None:
            method_versions.update(treeval_version.model_dump(exclude_unset=True, exclude_none=True))

    if method_versions:
        methods = SoftwareVersions(**method_versions)
    else:
        methods = {}

    # 3. Combine into a single  context
    context = {"taxonomy_info": taxonomy_info, "assembly": assembly_context, "methods": methods}
    logger.info("Generated context:")
    logger.info(context)

    # 4. Use the context to render the templates
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(Path(outdir, "rendered.md"), "wt") as f:
        f.write(render_template(args.template, **context))


if __name__ == "__main__":
    main()
