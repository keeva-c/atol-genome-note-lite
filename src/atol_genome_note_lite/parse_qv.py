#!/usr/bin/env python3

import csv
import json

path_to_field_mapping = "dev/qv_to_fields.csv"
path_to_file = "v"
field_mapping_dict = {}
qv_values_dict = {}
output_dict = {}
path_to_output = "dev/assembly_qv_fields.json"

#save genome note field names to a dictionary
with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[1]]=line[0]

#create a dictionary of qv values for primary, alt and combined assemblies
with open(path_to_file, "rt") as f:
    qv_table = csv.reader(f, delimiter='\t')
    header_row = next(qv_table)
    qv_position = header_row.index("QV") #define the position of the QV values in each row
    for row in qv_table:
        key = row[0]
        value = row[qv_position]
        qv_values_dict[key]=value

#map values from qv dictionary to genome note field names
for metadata_field,qv_row_name in field_mapping_dict.items():
    output_dict[metadata_field] = qv_values_dict[qv_row_name]

#save output to json file
with open(path_to_output, "wt") as f:
    json.dump(output_dict, f)