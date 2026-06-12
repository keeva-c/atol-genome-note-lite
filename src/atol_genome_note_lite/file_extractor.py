#!/usr/bin/env python3

import argparse
import json
import logging
import re
from pathlib import Path

# configure logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler("file_extractor.log", encoding="utf-8", mode="w")
logger.addHandler(console_handler)
logger.addHandler(file_handler)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)
logger.setLevel("DEBUG")
console_handler.setLevel("INFO")

# add arguments
argument_parser = argparse.ArgumentParser(
    description="This tool extracts file paths for files used in the genome note lite from a text file containing the directory contents generated in the genome assembly pipeline.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
input_group = argument_parser.add_argument_group("Input")
output_group = argument_parser.add_argument_group("Output")
input_group.add_argument(
    "--file_dir",
    type=Path,
    help="a text file listing all file paths contained in the genome assembly pipeline ouptut directory"
)
input_group.add_argument(
    "--tolid",
    help="the ToLID for the specimen used to generate the genome assembly"
)
argument_parser.add_argument(
    "--hic",
    action="store_true",
    help="extracts files for assemblies generated with Hi-C data"
)
output_group.add_argument(
    "--output",
    type=Path,
    default=Path("results/found_files.json"),
    help="a JSON file containing the list of extracted paths with their file types"
)
args = argument_parser.parse_args()

# TODO: handle output directories and file paths for assemblies generated without hi-c

# set global variables
tolid = args.tolid
file_paths = {}

hic_patterns = {
    'hap_1_summary_stats': f"{tolid}\\..+\\.phased/scaffolding/asm_hap1_scaffolds_final\\.fa\\.assembly_summary",
    'hap_2_summary_stats': f"{tolid}\\..+\\.phased/scaffolding/asm_hap2_scaffolds_final\\.fa\\.assembly_summary",
    'hap_1_busco_stats': f"{tolid}\\..+\\.phased/scaffolding/busco\\..+/short_summary\\.specific\\..+\\.asm_hap1_scaffolds_final\\.fa\\.json", # specifically taking the phased assembly busco
    'hap_2_busco_stats': f"{tolid}\\..+\\.phased/scaffolding/busco\\..+/short_summary\\.specific\\..+\\.asm_hap2_scaffolds_final\\.fa\\.json",
    'kmer_stats': f"{tolid}\\..+\\.phased/scaffolding/merqury\\..*/asm\\..*\\.completeness\\.stats",
    'qv_stats': f"{tolid}\\..+\\.phased/scaffolding/merqury\\..*/asm\\..*\\.qv",
    'hap_1_contact_map': f"{tolid}\\..+\\.phased/scaffolding/contact_maps/asm_hap1\\.pretext\\.FullMap\\.png",
    'hap_2_contact_map': f"{tolid}\\..+\\.phased/scaffolding/contact_maps/asm_hap2\\.pretext\\.FullMap\\.png",
    'mitogenome_stats': f"{tolid}\\.mitohifi/contigs_stats\\.tsv", # defaulting to mitohifi results generated from sequence reads (not assembled contigs)
    'software_versions': f"pipeline_info/genomeassembly_software_versions\\.yml",
    'genomescope_plot': f"kmer/k../long/{tolid}.long.k.._linear_plot\\.png" # currently missing in the results directory
}

no_hic_patterns = {
    'busco_stats': f"{tolid}\\..+\\.purged/purging/busco\\..+/short_summary\\.specific\\..+\\.purged\\.fa\\.json",
    'kmer_stats': f"{tolid}\\..+\\.purged/purging/merqury\\..*/asm\\..*\\.completeness\\.stats",
    'qv_stats': f"{tolid}\\..+\\.purged/purging/merqury\\..*/asm\\..*\\.qv",
    'summary_stats': f"{tolid}\\..+\\.purged/purging/asm\\.purged\\.fa\\.assembly_summary",
    'mitogenome_stats': f"{tolid}\\.mitohifi/contigs_stats\\.tsv", # defaulting to mitohifi results generated from sequence reads (not assembled contigs)
    'software_versions': f"pipeline_info/genomeassembly_software_versions\\.yml",
    'genomescope_plot': f"kmer/k../long/{tolid}.long.k.._linear_plot\\.png" # currently missing in the results directory
}

# functions
def find_file(file_type, pattern):
    logger.debug(f"Searching for {file_type}...")
    path_list = re.findall(pattern, file_names)
    if len(path_list) == 0:
        logger.warning(f"No file found for {file_type}.")
        path_list = None
    elif len(path_list) > 1:
        logger.warning(f"Multiple files found for {file_type}: {path_list}")
    else:
        logger.info(f"File path found for {file_type}")
        path_list = path_list[0]
#        path_list[0] = Path(path_list[0]) # convert string to Path
    return path_list

# configurations
if args.hic:
    patterns = hic_patterns
else:
    patterns = no_hic_patterns

# start script
with open(args.file_dir, "r", encoding="utf-8") as f:
    file_names = f.read()

for file_type, pattern in patterns.items():
    path = find_file(file_type, pattern)
    file_paths[file_type] = path

with open(args.output, "wt", encoding="utf-8") as f:
    logger.info(f"Writing file paths for {args.tolid} to {args.output}")
    json.dump(file_paths, f)

print(f"rclone copy command for {args.tolid}:")

for file_path in file_paths.values():
    if file_path:
        print(f"--include {file_path} `")