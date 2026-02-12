#!/usr/bin/env python3

import argparse
import jinja2
import json
import logging
import os
import sys
#from jinja2 import Template, Environment, FileSystemLoader, Undefined
from pathlib import Path

# configure logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler("genome_note_lite.log", encoding="utf-8", mode="w")
logger.addHandler(console_handler)
logger.addHandler(file_handler)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)
logger.setLevel("DEBUG")
console_handler.setLevel("INFO")

# setting the global default undefined value
class undefined_tokens(jinja2.Undefined):
    def __str__(self):
        return "*not provided*"

# set global variables
path_to_sample_supplement_output = "templates/supplementary_sample_data_for_genome_note.md"
path_to_extract_supplement_output = "templates/supplementary_extract_data_for_genome_note.md"
path_to_seq_supplement_output = "templates/supplementary_seq_data_for_genome_note.md"
path_to_bpa_package_supplement_output = "templates/supplementary_package_data_for_genome_note.md"
processed_wgs_file_paths = []
processed_hic_file_paths = []
input_sample_ids = []
input_library_ids = []

# setting the environment for the genome note templates
env = jinja2.Environment(loader=jinja2.FileSystemLoader("templates"),undefined=undefined_tokens)

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

# TODO: check that arguments are valid paths

# functions to preprocess metadata
def preprocess_metadata(metadata_file, processed_files):
    '''combining the above two metadata formatting functions together and saving output to a new file with '-processed' suffix'''
    processed_file = Path(metadata_file).stem + "-processed" + Path(metadata_file).suffix
    logger.info(f"Preprocessing {metadata_file} and saving ouput to interim {processed_file} file")
    with open(metadata_file, 'r') as f:
        unprocessed_file = json.load(f)
        metadata_w_initiative = map_bpa_initiative(unprocessed_file)
        metadata_w_platform = sanitise_platform(metadata_w_initiative)
        processed_metadata = overwrite_empty_strings(metadata_w_platform)
    with open(processed_file, 'w') as f:
        json.dump(processed_metadata, f)
    processed_files.append(processed_file)
    return processed_files

def map_bpa_initiative(metadata):
    '''add standardised bpa_initiative attribute based on bioplatforms_project_id metadata'''
    initiative_acronym = metadata.get('sample', {}).get('bioplatforms_project_id')
    if initiative_acronym == 'bpa-ipm':
        full_initiative = 'Integrated Pest Management Omics Initiative'
    elif initiative_acronym == 'threatened-species':
        full_initiative = 'Threatened Species Initiative'
    elif initiative_acronym == 'aus-fish':
        full_initiative = 'Australian Fish Genomics Initiative'
    elif initiative_acronym == 'aus-avian':
        full_initiative = 'Australian Avian Genomics Initiative'
    elif initiative_acronym == 'fungi':
        full_initiative = 'Australian Functional Fungi Initiative'
    elif initiative_acronym == 'cipps':
        full_initiative = 'ARC for Innovations in Peptide and Protein Science (CIPPS)'
    elif initiative_acronym == 'ausarg':
        full_initiative = 'Australian Amphibian and Reptile Genomics Initiative'
    elif initiative_acronym == 'forest-resilience':
        full_initiative = 'Genomics for Forest Resilience Initiative'
    elif initiative_acronym == 'bpa-plants':
        full_initiative = 'Genomics for Australian Plants'
    else:
        logger.warning(f"Bioplatforms initiative ({initiative_acronym}) not recognised.")
        full_initiative = None
    metadata['sample']['bpa_initiative'] = full_initiative
    return metadata

def sanitise_platform(metadata):
    '''sanitise the controlled vocab values required by INSDC to be more human-readable'''
    verbatim_platform = metadata.get('experiment', {}).get('platform')
    if verbatim_platform == 'PACBIO_SMRT':
        cleaned_platform = 'PacBio'
    elif verbatim_platform == 'OXFORD_NANOPORE':
        cleaned_platform = 'Oxford Nanopore'
    elif verbatim_platform == 'ILLUMINA':
        cleaned_platform = 'Illumina'
    else:
        logger.warning(f"Sequencing platform not recognised - platform will be rendered as {verbatim_platform}.")
        cleaned_platform = verbatim_platform
    metadata['experiment']['platform'] = cleaned_platform
    return(metadata)

def overwrite_empty_strings(metadata):
    '''updating metadata input to remove metadata containing empty strings'''
    output_dict = {}
    for database_sect, metadata_dict in metadata.items():
        updated_metadata = {}
        for field, value in metadata_dict.items():
            if value == "" or value is None:
                continue
            else:
                updated_metadata[field] = value
        output_dict[database_sect] = updated_metadata
    return output_dict

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

# functions to render helper files to append to core template
def render_if_not_duplicate(helper_file, helper_metadata, helper_sample, helper_library, enum_idx):
    '''checking if the sample and library ids in the metadata file have previously been encountered, and passing the metadata to the render_helper function if the sample or library has not perviously been rendered into markdown'''
    logger.info(f"Reading metadata from {helper_file} and rendering helper files")
    # render the sample helper file
    if helper_sample not in input_sample_ids:
        render_helper(
            helper_template="sup_sample_template.md",
            helper_metadata=helper_metadata,
            helper_output=path_to_sample_supplement_output,
            enum_idx=enum_idx
        )
    # render the extraction helper file
    if helper_library not in input_library_ids:
        render_helper(
            helper_template="sup_extract_template.md",
            helper_metadata=helper_metadata,
            helper_output=path_to_extract_supplement_output,
            enum_idx=enum_idx
        )
    # render the sequencing helper file
    render_helper(
        helper_template="sup_sequencing_template.md",
        helper_metadata=helper_metadata,
        helper_output=path_to_seq_supplement_output,
        enum_idx=enum_idx
    )
    # render the package helper file
    render_helper(
        helper_template="sup_bpa_package_template.md",
        helper_metadata=helper_metadata,
        helper_output=path_to_bpa_package_supplement_output,
        enum_idx=enum_idx
    )

def render_helper(helper_template, helper_metadata, helper_output, enum_idx):
    '''rendering the supplementary helper files, and appending the newly rendered markdown to the original markdown if multiple metadata files were provided in the original input arguments'''
    sup_template = env.get_template(helper_template)
    sup_render = sup_template.render(helper_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
    if enum_idx == 1 or not os.path.exists(helper_output):
        with open(helper_output, "wt", encoding="utf-8") as f:
            f.write(sup_render)
    elif enum_idx > 1:
        with open(helper_output, "at", encoding="utf-8") as f:
            if helper_template != "sup_bpa_package_template.md":
                f.write('\n' + sup_render)
            elif helper_template == "sup_bpa_package_template.md":
                f.write(sup_render)

logger.info("Starting script")

# preprocessing metadata for input WGS metadata
for input_file in args.wgs_metadata:
    preprocess_metadata(input_file, processed_wgs_file_paths) 

# preprocessing metadata for input Hi-C 
if args.hic_metadata is not None:
    for input_file in args.hic_metadata:
        preprocess_metadata(input_file, processed_hic_file_paths) 

# initialise lists for cross-checking sample and library ids (to avoid duplicating metadata in the genome note lite)
all_input_files = list(processed_wgs_file_paths)

if args.hic_metadata is not None:
    all_input_files.extend(processed_hic_file_paths)

# read Hi-C and supplementary WGS input metadata and render supplementary helper files
for idx, file in enumerate(all_input_files):
    # find the sample and library IDs and append to a list
    with open(file, "rt") as f:
        full_metadata = json.load(f)
        sample_name = full_metadata['sample']['bpa_sample_id']
        library_name = full_metadata['experiment']['bpa_library_id']
        if idx == 0:
            next(enumerate(all_input_files))
        else:
            render_if_not_duplicate(file, full_metadata, sample_name, library_name, idx)
    logger.debug(f"Appending {sample_name} to input sample IDs list from {file}")
    input_sample_ids.append(sample_name)
    logger.debug(f"Appending {library_name} to input library IDs list from {file}")
    input_library_ids.append(library_name)

# initialise main template
if args.w_annotation:
    logger.info("Preparing genome note lite template for annotated assembly")
    template = env.get_template("genome-note-lite-annot-template.md")
elif args.wo_annotation:
    logger.info("Preparing genome note lite template for assembly without annotation")
    template = env.get_template("genome-note-lite-asm-only-template.md")
else:
    logger.info("Preparing genome note lite template for assembly without annotation (by default)")
    template = env.get_template("genome-note-lite-asm-only-template.md")

# prepare output directory
if str(args.output) == "results/genome_note_lite.md":
    if not os.path.isdir("results"):
        os.mkdir("results")

# read the sample metadata
logger.info(f"Reading metadata from {processed_wgs_file_paths[0]}")
with open(processed_wgs_file_paths[0], "rt") as f:
    assembly_sample_metadata = json.load(f)

# render the core template, integrating the supplementary markdown files for hi-c and/or secondary wgs data
logger.info(f"Rendering core template, including helper files (if applicable) and writing output to {args.output}")
with open(args.output, "wt", encoding="utf-8") as f:
    try:
        f.write(template.render(assembly_sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal,add_spaces=add_spaces))
    except jinja2.exceptions.UndefinedError as e_jinja_undefined:
        logger.exception(f"Some expected metadata object (e.g. organism, sample, experiment, runs, assembly, or annotation) was not found in the input metadata: {e_jinja_undefined}")

# remove processed metadata and supplementary markdown files
logger.info("Removing processed metadata and supplementary helper files and ending script")
for file in all_input_files:
    os.remove(file)

for file in [path_to_sample_supplement_output, path_to_extract_supplement_output, path_to_seq_supplement_output, path_to_bpa_package_supplement_output]:
    if os.path.exists(file):
        os.remove(file)