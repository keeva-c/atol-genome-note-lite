#!/usr/bin/env python3

import csv
import json

path_to_field_mapping = "dev/kmer_to_fields.csv"
path_to_file = " "
field_mapping_dict = {}
kmer_values_dict = {}
output_dict = {}

#save genome note field names to a dictionary
with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[1]]=line[0]

#create a dictionary of kmer completeness values for primary, alt and combined assemblies
with open(path_to_file, "rt") as f:
    kmer_table = csv.reader(f, delimiter='\t')
    header_row = next(kmer_table)
    kmer_position = header_row.index("% Covered") #define the position of the kmer completeness values in each row
    for row in kmer_table:
        key = row[0]
        value = row[kmer_position]
        kmer_values_dict[key]=value

#map values from kmer completeness dictionary to genome note field names
for metadata_field,kmer_row_name in field_mapping_dict.items():
    output_dict[metadata_field] = kmer_values_dict[kmer_row_name]

print(json.dumps(output_dict))