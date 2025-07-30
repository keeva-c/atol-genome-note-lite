#!/usr/bin/env python3

import csv
import json

path_to_field_mapping = "dev/read_stats_to_fields.csv"
path_to_file = ""
field_mapping_dict = {}
SN_initial = "SN	"
SN_line_list = []
read_stat_dict = {}
output_dict = {}
path_to_output = "dev/run_stats_fields.json"

#save genome note field names to a dictionary
with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[1]]=line[0]

with open(path_to_file, "rt") as f:
    for line in f:
        if SN_initial in line:
            SN_lines = line.strip().split('\t')
            SN_line_list.append(SN_lines) #create a list of lines containing summary number stats
            for SN_line in SN_line_list:
                key = SN_line[1]
                value = SN_line[2]
                read_stat_dict[key] = value #create  a dictionary of read stat field names and values

#map values from read stats dictionary to genome note field names
for metadata_field,read_stat_field in field_mapping_dict.items():
    output_dict[metadata_field] = read_stat_dict[read_stat_field]

#print json output
print(json.dumps(output_dict))

#write output to json file
with open(path_to_output, "wt") as f:
    json.dump(output_dict, f)