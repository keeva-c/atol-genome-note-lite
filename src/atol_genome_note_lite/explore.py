#!/usr/bin/env python3

import json
import argparse
from pathlib import Path
from jinja2 import Template, Environment, FileSystemLoader, Undefined

print("Starting script")

#path_to_sample_metadata = "dev/sample_WGS.json"
#path_to_WGS_supplement_metadata = None #set to None if n/a
#path_to_hic_supplement_metadata = "dev/sample_hic.json" #set to None if n/a
path_to_sample_supplement_output = "dev/supplementary_sample_data_for_genome_note.md"
path_to_extract_supplement_output = "dev/supplementary_extract_data_for_genome_note.md"
path_to_seq_supplement_output = "dev/supplementary_seq_data_for_genome_note.md"
path_to_bpa_package_supplement_output = "dev/supplementary_package_data_for_genome_note.md"
#path_to_genome_note_output = "results/genome_note_lite.md"

argument_parser = argparse.ArgumentParser(
    description="This tool generates a draft genome note lite markdown document based on sample, read, assembly, and annotation metadata."
)
input_group = argument_parser.add_argument_group("Input")
output_group = argument_parser.add_argument_group("Output")
input_group.add_argument(
    "--wgs_metadata",
    nargs="+",
    type=Path,
    help="at least one JSON file listing metadata for a PacBio WGS sample and sequencing run. The first file should also contain assembly metadata and optionally annotation metadata."
)
input_group.add_argument(
    "--hic_metadata",
    nargs="?",
    type=Path,
    help="an optional JSON file listing metadata for a Hi-C sample and sequencing run."
)
output_group.add_argument(
    "--output",
    type=Path,
    default=Path("results/genome_note_lite.md"),
    help="optional output file address."
)
argument_parser.add_argument(
    "--w_annotation",
    action="store_true",
    help="runs the genome note lite for an assembly with an annotation."
)
argument_parser.add_argument(
    "--wo_annotation",
    action="store_true",
    help="runs the genome note lite for an assembly only (no annotation)."
)
args = argument_parser.parse_args()

# updating metadata input to remove metadata containing empty strings
def overwrite_empty_strings(metadata):
    output_dict = {}
    for dict_name, key_value in metadata.items():
        updated_metadata = {}
        for key, value in key_value.items():
            if value == "" or value is None:
                updated_metadata = updated_metadata
            else:
                updated_metadata[key] = value
        output_dict[dict_name] = updated_metadata
    return output_dict

# add standardised bpa_initiative attribute based on bioplatforms_project_id metadata
def map_bpa_initiative(metadata):
    for dict_name, key_value in metadata.items():
        if dict_name == 'sample':
            if key_value['bioplatforms_project_id'] == 'bpa-ipm':
                key_value['bpa_initiative'] = 'Integrated Pest Management Omics Initiative'
            elif key_value['bioplatforms_project_id'] == 'threatened-species':
                key_value['bpa_initiative'] = 'Threatened Species Initiative'
            elif key_value['bioplatforms_project_id'] == 'aus-fish':
                key_value['bpa_initiative'] = 'Australian Fish Genomics Initiative'
            elif key_value['bioplatforms_project_id'] == 'plant-pathogen':
                key_value['bpa_initiative'] = 'Plant Pathogen Omics Initiative'
            elif key_value['bioplatforms_project_id'] == 'aus-avian':
                key_value['bpa_initiative'] = 'Australian Avian Genomics Initiative'
            elif key_value['bioplatforms_project_id'] == 'fungi':
                key_value['bpa_initiative'] = 'Australian Functional Fungi Initiative'
            elif key_value['bioplatforms_project_id'] == 'cipps':
                key_value['bpa_initiative'] = 'ARC for Innovations in Peptide and Protein Science (CIPPS)'
            elif key_value['bioplatforms_project_id'] == 'ausarg':
                key_value['bpa_initiative'] = 'Australian Amphibian and Reptile Genomics Initiative'
            elif key_value['bioplatforms_project_id'] == 'forest-resilience':
                key_value['bpa_initiative'] = 'Genomics for Forest Resilience Initiative'
            elif key_value['bioplatforms_project_id'] == 'grasslands':
                key_value['bpa_initiative'] = 'Australian Grasslands Initiative'
            elif key_value['bioplatforms_project_id'] == 'bpa-plants':
                key_value['bpa_initiative'] = 'Genomics for Australian Plants'
            else:
                key_value['bpa_initiative'] = None
    return metadata

print("Preprocessing metadata input")

# preprocessing metadata for input WGS metadata
with open(args.wgs_metadata[0], 'r') as f:
    unprocessed_sample_metadata = json.load(f)

sample_metadata_w_initiative = map_bpa_initiative(unprocessed_sample_metadata)
processed_sample_metadata = overwrite_empty_strings(sample_metadata_w_initiative)

with open(args.wgs_metadata[0], 'w') as f:
    json.dump(processed_sample_metadata, f)

# preprocessing metadata for secondary input WGS metadata file if provided
if len(args.wgs_metadata) > 1:
    with open(args.wgs_metadata[1], 'r') as f:
        unprocessed_wgs_sup_metadata = json.load(f)

    wgs_sup_metadata_w_initiative = map_bpa_initiative(unprocessed_wgs_sup_metadata)
    processed_wgs_sup_metadata = overwrite_empty_strings(wgs_sup_metadata_w_initiative)

    with open(args.wgs_metadata[1], 'w') as f:
        json.dump(processed_wgs_sup_metadata, f)

# preprocessing metadata for input Hi-C metadata if provided
if args.hic_metadata is not None:
    with open(args.hic_metadata, 'r') as f:
        unprocessed_hic_sup_metadata = json.load(f)

    hic_sup_metadata_w_initiative = map_bpa_initiative(unprocessed_hic_sup_metadata)
    processed_hic_sup_metadata = overwrite_empty_strings(hic_sup_metadata_w_initiative)

    with open(args.hic_metadata, 'w') as f:
        json.dump(processed_hic_sup_metadata, f)

# setting the global default undefined value
class undefined_tokens(Undefined):
    def __str__(self):
        return "*not provided*"

# functions to make things in the template pretty
def make_pretty_number(ugly_number):
    if float(ugly_number) > 1_000:
        return "{:,}".format(int(ugly_number))
    else:
        return ugly_number

def round_decimal(decimal):
    if decimal is not None:
        return "{:,.1f}".format(float(decimal))
    else:
        return decimal

def round_bases_up(base_pairs):
    if int(base_pairs) > 1_000_000_000:
        return "{:.2f} Gb".format(int(base_pairs) / 1_000_000_000)
    elif int(base_pairs) > 1_000_000:
        return "{:.2f} Mb".format(int(base_pairs) / 1_000_000)
    elif int(base_pairs) > 1_000:
        return "{:.2f} kb".format(int(base_pairs) / 1_000)
    else:
        return base_pairs + " bp"

def add_spaces(long_string):
    separator = " "
    formatted_string = ""
    for character in long_string:
        if character != ",":
            formatted_string = formatted_string + character
        elif character == ",":
            formatted_string = formatted_string + character + separator
    return formatted_string

# setting the environment for the genome note templates
env = Environment(loader=FileSystemLoader("templates"),undefined=undefined_tokens)

# check whether hi-c metadata are provided
if args.hic_metadata is not None:

    # find the sample IDs for WGS and Hi-C packages
    with open(args.wgs_metadata[0], "rt") as f:
        sample_wgs_metadata = json.load(f)
        sample_wgs_name = sample_wgs_metadata['sample']['bpa_sample_id']

    with open(args.hic_metadata, "rt") as f:
        sample_hic_metadata = json.load(f)
        sample_hic_name = sample_hic_metadata['sample']['bpa_sample_id']

    # read the supplementary templates
    # only read supplementary sample templates if the WGS and Hi-C packages are derived from different samples
    if sample_wgs_name != sample_hic_name:
        sup_samp_template = env.get_template("sup_sample_template.md")

    sup_ext_template = env.get_template("sup_extract_template.md")
    sup_seq_template = env.get_template("sup_sequencing_template.md")
    sup_package_template = env.get_template("sup_bpa_package_template.md")

    # read the supplementary sample metadata
    print(f"Reading Hi-C supplementary sample metadata from {args.hic_metadata}")
    with open(args.hic_metadata, "rt") as f:
        supplement_metadata = json.load(f)

    # render the supplementary templates
    print(f"Writings Hi-C supplementary output to helper files")
    # only render the supplementary sample templates if the WGS and Hi-C packages are derived from different samples
    if sample_wgs_name != sample_hic_name:
        sup_sample_render = sup_samp_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
        with open(path_to_sample_supplement_output, "wt", encoding="utf-8") as f:
            f.write(sup_sample_render)

    sup_ext_render = sup_ext_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
    sup_seq_render = sup_seq_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
    sup_package_render = sup_package_template.render(supplement_metadata)

    with open(path_to_extract_supplement_output, "wt", encoding="utf-8") as f:
        f.write(sup_ext_render)
    with open(path_to_seq_supplement_output, "wt", encoding="utf-8") as f:
        f.write(sup_seq_render)
    with open(path_to_bpa_package_supplement_output, "wt", encoding="utf-8") as f:
        f.write(sup_package_render)

else:
    print("No Hi-C data provided; genome note lite script running with WGS data only")

# check whether supplementary wgs metadata are provided
if len(args.wgs_metadata) > 1:

    # find the sample and library IDs for WGS and supplementary WGS packages
    with open(args.wgs_metadata[0], "rt") as f:
        sample_1_metadata = json.load(f)
        sample_1_name = sample_1_metadata['sample']['bpa_sample_id']
        sample_1_library = sample_1_metadata['experiment']['bpa_library_id']

    with open(args.wgs_metadata[1], "rt") as f:
        sample_2_metadata = json.load(f)
        sample_2_name = sample_2_metadata['sample']['bpa_sample_id']
        sample_2_library = sample_2_metadata['experiment']['bpa_library_id']

    # read the supplementary templates
    # only read supplementary sample templates if the WGS packages are derived from different samples
    if sample_1_name != sample_2_name:
        sup_samp_template = env.get_template("sup_sample_template.md")

    # only read supplementary extract templates if the WGS packages are derived from different libraries
    if sample_1_library != sample_2_library:
        sup_ext_template = env.get_template("sup_extract_template.md")

    sup_seq_template = env.get_template("sup_sequencing_template.md")
    sup_package_template = env.get_template("sup_bpa_package_template.md")

    # read the supplementary sample metadata
    print(f"Reading WGS supplementary sample metadata from {args.wgs_metadata[1]}")
    with open(args.wgs_metadata[1], "rt") as f:
        WGS_supplement_metadata = json.load(f)

    # render the supplementary templates
    print(f"Writings WGS supplementary output to helper files")
    # only render the supplementary sample templates if the WGS packages are derived from different samples
    if sample_1_name != sample_2_name:
        sup_sample_render = sup_samp_template.render(WGS_supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
        with open(path_to_sample_supplement_output, "rt", encoding="utf-8") as f:
            sample_supplement_md = f.read()
        combined_sample_sup = sup_sample_render  + '\n' + sample_supplement_md
        with open(path_to_sample_supplement_output, "wt", encoding="utf-8") as f:
            f.write(combined_sample_sup)

    # only render the supplementary extraction templates if the WGS packages are derived from different libraries
    if sample_1_library != sample_2_library:
        sup_ext_render = sup_ext_template.render(WGS_supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
        with open(path_to_extract_supplement_output, "rt", encoding="utf-8") as f:
            extract_supplement_md = f.read()
        combined_extract_sup = sup_ext_render + '\n' + extract_supplement_md
        with open(path_to_extract_supplement_output, "wt", encoding="utf-8") as f:    
            f.write(combined_extract_sup)

    sup_seq_render = sup_seq_template.render(WGS_supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
    with open(path_to_seq_supplement_output, "rt", encoding="utf-8") as f:
        sequence_supplement_md = f.read()
    combined_sequence_sup = sup_seq_render + '\n' + sequence_supplement_md
    with open(path_to_seq_supplement_output, "wt", encoding="utf-8") as f:
        f.write(combined_sequence_sup)

    sup_package_render = sup_package_template.render(WGS_supplement_metadata)
    with open(path_to_bpa_package_supplement_output, "rt", encoding="utf-8") as f:
        package_supplement_md = f.read()
    combined_package_sup = sup_package_render + package_supplement_md
    with open(path_to_bpa_package_supplement_output, "wt", encoding="utf-8") as f:
        f.write(combined_package_sup)

else:
    print("One WGS package provided")

# initialise template
if args.w_annotation:
    print("Preparing genome note lite for annotated assembly")
    template = env.get_template("genome-note-lite-annot-template.md")
elif args.wo_annotation:
    print("Preparing genome note lite for assembly without annotation")
    template = env.get_template("not-a-genome-note-lite-template.md")
else:
    print("Preparing genome note lite for assembly without annotation (by default)")
    template = env.get_template("not-a-genome-note-lite-template.md")

# read the sample metadata
print(f"Reading sample metadata from {args.wgs_metadata[0]}")
with open(args.wgs_metadata[0], "rt") as f:
    sample_metadata = json.load(f)

# render the core template, integrating the supplementary markdown files for hi-c and/or secondary wgs data
print(f"Combining and writing output to {args.output}")
with open(args.output, "wt", encoding="utf-8") as f:
    f.write(template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal,add_spaces=add_spaces))

#wipe supplementary markdown files
print("Wiping supplementary helper files and ending script")
with open(path_to_sample_supplement_output, "wt") as f:
    f.close()
with open(path_to_extract_supplement_output, "wt") as f:
    f.close()
with open(path_to_seq_supplement_output, "wt") as f:
    f.close()
with open(path_to_bpa_package_supplement_output, "wt") as f:
    f.close()