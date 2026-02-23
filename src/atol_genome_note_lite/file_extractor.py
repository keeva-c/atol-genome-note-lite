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
argument_parser.add_argument(
    "--file_dir",
    type=Path,
    help="A text file listing all file paths contained in the genome assembly pipeline ouptut directory."
)
argument_parser.add_argument(
    "--tolid",
    help="The ToLID for the specimen used to generate the genome assembly."
)
argument_parser.add_argument(
    "--output",
    type=Path,
    default=Path("results/found_files.json"),
    help="A JSON file containing the list of extracted paths with their file types."
)
args = argument_parser.parse_args()

# TODO: handle output directories and file paths for assemblies generated without hi-c

# set global variables
tolid = args.tolid
file_paths = {}

patterns = {
    'busco_stats': f"{tolid}\\.hifiasm-hic.*/scaffolding_hap1/yahs/asm_hap1_scaffolds_final.*short_summary.json",
    'kmer_stats': f"{tolid}\\.hifiasm-hic.*/scaffolding_hap1/yahs/asm_hap1_scaffolds_final\\.fa\\.ccs\\.merquryfk/{tolid}\\.ccs\\.completeness\\.stats",
    'qv_stats': f"{tolid}\\.hifiasm-hic.*/scaffolding_hap1/yahs/asm_hap1_scaffolds_final\\.fa\\.ccs\\.merquryfk/{tolid}\\.ccs\\.qv",
    'summary_stats': f"{tolid}\\.hifiasm-hic.*/scaffolding_hap1/yahs/asm_hap1_scaffolds_final\\.fa\\.gz\\.assembly_summary",
    'contact_map': f"{tolid}\\.hifiasm-hic.*/scaffolding_hap1/yahs/asm_hap1\\.pretext\\.FullMap\\.png",
    'mitogenome_stats': f"{tolid}\\.hifiasm-hic.*/mito/contigs_stats\\.tsv",
    'software_versions': f"pipeline_info/genomeassembly_software_versions\\.yml",
    'genomescope_plot': f"kmer/k../long/{tolid}.long.k.._linear_plot\\.png"
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

# start script
with open(args.file_dir, "r") as f:
    file_names = f.read()

for file_type, pattern in patterns.items():
    path = find_file(file_type, pattern)
    file_paths[file_type] = path

with open(args.output, "wt", encoding="utf-8") as f:
    logger.info(f"Writing file paths for {args.tolid} to {args.output}")
    json.dump(file_paths, f)