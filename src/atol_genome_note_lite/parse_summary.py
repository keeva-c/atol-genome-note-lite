#!/usr/bin/env python3

import csv
import json

path_to_field_mapping = "dev/summary_to_fields.csv"
path_to_file = ""
sep = ": "
values_dict = {}
output_dict = {}
field_mapping_dict = {}

path_to_output = "dev/assembly_summary_fields.json"

with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[1]]=line[0]

with open(path_to_file, "rt") as f:
    next(f) #take out the header
    for line in f:
        splits = line.strip().split(sep)
        if len(splits) == 2:
            key = splits[0]
            value = splits[1]
            values_dict[key]=value

for metadata_field,summary_field_name in field_mapping_dict.items():
    output_dict[metadata_field] = values_dict[summary_field_name]

#write output to json file
with open(path_to_output, "wt") as f:
    json.dump(output_dict, f)