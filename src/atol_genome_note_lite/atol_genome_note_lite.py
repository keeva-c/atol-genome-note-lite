#!/usr/bin/env python3

import argparse
from importlib.resources import files
import json
from pathlib import Path

import jinja2

from atol_genome_note_lite.common import _existing_file, logger
from atol_genome_note_lite.models import TaxonomyInfo
from canopy_client import CanopySession
from yaml_manifest import Manifest


class IgnoreUndefined(jinja2.Undefined):
    def _fail_with_undefined_error(self, *args, **kwargs):
        logger.warning(self._undefined_message)
        return None


def parse_arguments() -> argparse.Namespace:
    my_files = files()

    parser = argparse.ArgumentParser(
        description="Render the templates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    _ = parser.add_argument(
        "--template",
        type=_existing_file,
        default=Path(my_files, "templates", "base.md.j2"),
        help="Base template.",
    )

    _ = parser.add_argument(
        "--outdir", type=Path, help="Output directory.", required=True
    )

    _ = parser.add_argument(
        "manifest", type=_existing_file, help="AToL manifest for this assembly."
    )

    return parser.parse_args()


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


def seqstats():
    pass


def main():

    args = parse_arguments()
    logger.name = "atol-genome-note-lite"

    # 0. Read the assembly Manifest
    with open(args.manifest, "rb") as f:
        manifest = Manifest.model_validate_json(f.read())
    taxon_id = manifest.taxon_id

    # 1. Read the required information off Canopy
    canopy_session = CanopySession()

    canopy_taxonomy_info = canopy_session.get_taxonomy_info(taxon_id=taxon_id).json()
    taxonomy_info = TaxonomyInfo(**canopy_taxonomy_info)

    # 2. Harvest info from results files

    # 3. Combine into a single  context
    context = {"taxonomy_info": taxonomy_info.model_dump()}

    # 4. Use the context to render the templates
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(Path(outdir, "rendered.md"), "wt") as f:
        f.write(render_template(args.template, **context))


if __name__ == "__main__":
    main()
