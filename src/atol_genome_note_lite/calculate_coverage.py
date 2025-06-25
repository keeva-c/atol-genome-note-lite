#!/usr/bin/env python3

import json

path_to_read_stats = "dev/run_stats_fields.json"
path_to_assembly_stats = "dev/assembly_summary_fields.json"
coverage_dict = {}
path_to_output = "dev/assembly_coverage_field.json"

# retrieve base count
with open(path_to_read_stats, "rt") as f:
    read_stats = json.load(f)
    base_count = read_stats['run_base_count']

# retrieve genome size
with open(path_to_assembly_stats, "rt") as f:
    assembly_stats = json.load(f)
    genome_size = assembly_stats['genome_length']

# calculate coverage
coverage = int(base_count) / int(genome_size)

# save coverage as dictionary
coverage_dict["coverage"] = coverage

#save dictionary to json file
with open(path_to_output, "wt") as f:
    json.dump(coverage_dict, f)