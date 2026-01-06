#!/usr/bin/env python3

import argparse
import json
import os
from jinja2 import Template, Environment, FileSystemLoader, Undefined
from pathlib import Path

# set global variables
path_to_sample_supplement_output = "templates/supplementary_sample_data_for_genome_note.md"
path_to_extract_supplement_output = "templates/supplementary_extract_data_for_genome_note.md"
path_to_seq_supplement_output = "templates/supplementary_seq_data_for_genome_note.md"
path_to_bpa_package_supplement_output = "templates/supplementary_package_data_for_genome_note.md"

# set input arguments
argument_parser = argparse.ArgumentParser(
    description="This tool generates a draft genome note lite markdown document based on sample, read, assembly, and annotation metadata.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
input_group = argument_parser.add_argument_group("Input")
output_group = argument_parser.add_argument_group("Output")
input_group.add_argument(
    "--wgs_metadata",
    nargs="+",
    type=Path,
    help="at least one JSON file listing metadata for a WGS sample and sequencing run. The first file should also contain assembly metadata and optionally annotation metadata."
)
input_group.add_argument(
    "--hic_metadata",
    nargs="*",
    type=Path,
    help="optional JSON file/s listing metadata for a Hi-C sample and sequencing run."
)
output_group.add_argument(
    "--output",
    default=Path("results/genome_note_lite.md"),
    type=Path,
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

print("Starting script")

# TODO: check that arguments are valid paths

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

processed_wgs_file_paths = []
processed_hic_file_paths = []

# preprocessing metadata for input WGS metadata
for input_file in args.wgs_metadata:
    with open(input_file, 'r') as f:
        unprocessed_sample_metadata = json.load(f)
        sample_metadata_w_initiative = map_bpa_initiative(unprocessed_sample_metadata)
        processed_sample_metadata = overwrite_empty_strings(sample_metadata_w_initiative)
    processed_file = Path(input_file).stem + "-processed" + Path(input_file).suffix

    with open((processed_file), 'w') as f:
        json.dump(processed_sample_metadata, f)
    processed_wgs_file_paths.append(processed_file)

# preprocessing metadata for input Hi-C metadata if provided
if args.hic_metadata is not None:
    for input_file in args.hic_metadata:
        with open(input_file, 'r') as f:
            unprocessed_hic_sup_metadata = json.load(f)
            hic_sup_metadata_w_initiative = map_bpa_initiative(unprocessed_hic_sup_metadata)
            processed_hic_sup_metadata = overwrite_empty_strings(hic_sup_metadata_w_initiative)
        processed_file = Path(input_file).stem + "-processed" + Path(input_file).suffix

        with open((processed_file), 'w') as f:
            json.dump(processed_hic_sup_metadata, f)
        processed_hic_file_paths.append(processed_file)

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

# setting the global default undefined value
class undefined_tokens(Undefined):
    def __str__(self):
        return "*not provided*"

# setting the environment for the genome note templates
env = Environment(loader=FileSystemLoader("templates"),undefined=undefined_tokens)

# initialise lists for cross-checking sample and library ids (to avoid duplicating metadata in the genome note lite)
all_input_files = list(processed_wgs_file_paths)

if args.hic_metadata is not None:
    all_input_files.extend(processed_hic_file_paths)

input_sample_ids = []
input_library_ids = []

# read Hi-C and supplementary WGS input metadata and render supplementary helper files
for idx, file in enumerate(all_input_files):
    # find the sample IDs and append to a list
    with open(file, "rt") as f:
        sample_metadata = json.load(f)
        sample_name = sample_metadata['sample']['bpa_sample_id']
        library_name = sample_metadata['experiment']['bpa_library_id']
        if idx == 1:
            print(f"Reading metadata from", file, "and rendering helper files")
            # render the sample helper file
            if sample_name not in input_sample_ids:
                sup_samp_template = env.get_template("sup_sample_template.md")
                sup_sample_render = sup_samp_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
                with open(path_to_sample_supplement_output, "wt", encoding="utf-8") as f:
                    f.write(sup_sample_render)
            # render the extraction helper file
            if library_name not in input_library_ids:
                sup_ext_template = env.get_template("sup_extract_template.md")
                sup_ext_render = sup_ext_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
                with open(path_to_extract_supplement_output, "wt", encoding="utf-8") as f:
                    f.write(sup_ext_render)
            # render the sequencing helper file
            sup_seq_template = env.get_template("sup_sequencing_template.md")
            sup_seq_render = sup_seq_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
            with open(path_to_seq_supplement_output, "wt", encoding="utf-8") as f:
                f.write(sup_seq_render)
            # render the package helper file
            sup_package_template = env.get_template("sup_bpa_package_template.md")
            sup_package_render = sup_package_template.render(sample_metadata)
            with open(path_to_bpa_package_supplement_output, "wt", encoding="utf-8") as f:
                f.write(sup_package_render)
        elif idx > 1:
            print(f"Reading metadata from", file, "and rendering helper files")
            # render the sample helper file
            if sample_name not in input_sample_ids:
                sup_samp_template = env.get_template("sup_sample_template.md")
                sup_sample_render = sup_samp_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
                with open(path_to_sample_supplement_output, "rt", encoding="utf-8") as f:
                    sample_supplement_md = f.read()
                combined_sample_sup = sample_supplement_md + '\n' + sup_sample_render
                with open(path_to_sample_supplement_output, "wt", encoding="utf-8") as f:
                    f.write(combined_sample_sup)
            # render the extraction helper file
            if library_name not in input_library_ids:
                sup_ext_template = env.get_template("sup_extract_template.md")
                sup_ext_render = sup_ext_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
                with open(path_to_extract_supplement_output, "rt", encoding="utf-8") as f:
                    extract_supplement_md = f.read()
                combined_extract_sup = extract_supplement_md + '\n' + sup_ext_render
                with open(path_to_extract_supplement_output, "wt", encoding="utf-8") as f:
                    f.write(combined_extract_sup)
            # render the sequencing helper file
            sup_seq_template = env.get_template("sup_sequencing_template.md")
            sup_seq_render = sup_seq_template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
            with open(path_to_seq_supplement_output, "rt", encoding="utf-8") as f:
                sequence_supplement_md = f.read()
            combined_sequence_sup = sequence_supplement_md + '\n' + sup_seq_render
            with open(path_to_seq_supplement_output, "wt", encoding="utf-8") as f:
                f.write(combined_sequence_sup)
            # render the package helper file
            sup_package_template = env.get_template("sup_bpa_package_template.md")
            sup_package_render = sup_package_template.render(sample_metadata)
            with open(path_to_bpa_package_supplement_output, "rt", encoding="utf-8") as f:
                package_supplement_md = f.read()
            combined_package_sup = package_supplement_md + sup_package_render
            with open(path_to_bpa_package_supplement_output, "wt", encoding="utf-8") as f:
                f.write(combined_package_sup)
        else:
            next(enumerate(all_input_files))
    input_sample_ids.append(sample_name)
    input_library_ids.append(library_name)

# initialise main template
if args.w_annotation:
    print("Preparing genome note lite for annotated assembly")
    template = env.get_template("genome-note-lite-annot-template.md")
elif args.wo_annotation:
    print("Preparing genome note lite for assembly without annotation")
    template = env.get_template("not-a-genome-note-lite-template.md")
else:
    print("Preparing genome note lite for assembly without annotation (by default)")
    template = env.get_template("not-a-genome-note-lite-template.md")

# prepare output directory
if str(args.output) == "results/genome_note_lite.md":
    if not os.path.isdir("results"):
        os.mkdir("results")

# read the sample metadata
print(f"Reading metadata from {processed_wgs_file_paths[0]}")
with open(processed_wgs_file_paths[0], "rt") as f:
    assembly_sample_metadata = json.load(f)

# render the core template, integrating the supplementary markdown files for hi-c and/or secondary wgs data
print(f"Combining and writing output to {args.output}")
with open(args.output, "wt", encoding="utf-8") as f:
    f.write(template.render(assembly_sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal,add_spaces=add_spaces))

# remove processed metadata and supplementary markdown files
print("Removing processed metadata and supplementary helper files and ending script")
for file in all_input_files:
    os.remove(file)

for file in [path_to_sample_supplement_output, path_to_extract_supplement_output, path_to_seq_supplement_output, path_to_bpa_package_supplement_output]:
    if os.path.exists(file):
        os.remove(file)