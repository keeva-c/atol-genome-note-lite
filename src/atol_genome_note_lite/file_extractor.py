#!/usr/bin/env python3

import json
import logging
import re
# from pathlib import Path

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

# set global variables
file_dir = # file directory path
tolid = # tolid
# hic = True
file_paths = {}
output_json = "results/found_files.json"

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

def find_file(file_type, pattern):
    logger.debug(f"Searching for {file_type}...")
    path_list = re.findall(pattern, file_names)
    if len(path_list) == 0:
        logger.warning(f"No file found for {file_type}.")
        path_list.append(None)
    elif len(path_list) > 1:
        logger.warning(f"Multiple files found for {file_type}: {path_list}")
    else:
        logger.info(f"File path found for {file_type}")
        path_list = path_list[0]
#        path_list[0] = Path(path_list[0]) # convert string to Path
    return path_list

with open(file_dir, "r") as f:
    file_names = f.read()

for file_type, pattern in patterns.items():
    path = find_file(file_type, pattern)
    file_paths[file_type] = path

with open(output_json, "wt", encoding="utf-8") as f:
    logger.info(f"Writing output to {output_json}")
    json.dump(file_paths, f)