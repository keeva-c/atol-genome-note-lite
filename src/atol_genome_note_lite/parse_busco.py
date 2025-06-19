#!/usr/bin/env python3

import json
import csv

path_to_field_mapping = "dev/busco_to_fields.csv"
path_to_file = " "
field_mapping_dict = {}
busco_results = {}
output_dict = {}

#save genome note field names to a dictionary
with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[1]]=line[0]

#extract "results" and "lineage_dataset" objects from json file and merge into one dictionary
with open(path_to_file, "rt") as f:
    busco_stats = json.load(f)
    busco_results = busco_stats['results']
    busco_lineage = busco_stats['lineage_dataset']
    busco_for_mapping = busco_results | busco_lineage

#map values from merged dictionary to genome note field names
for metadata_field,busco_field_name in field_mapping_dict.items():
    output_dict[metadata_field] = busco_for_mapping[busco_field_name]

print(json.dumps(output_dict))