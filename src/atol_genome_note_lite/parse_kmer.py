#!/usr/bin/env python3

import csv
import json

path_to_field_mapping = "dev/kmer_to_fields.csv"
path_to_file = ""
field_mapping_dict = {}
kmer_values_dict = {}
output_dict = {}
path_to_output = "dev/assembly_kmer_fields.json"

#save genome note field names to a dictionary
with open(path_to_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        field_mapping_dict[line[0]]=line[1]

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
for kmer_row_index,metadata_field in field_mapping_dict.items():
    list_of_kmer_values = list(kmer_values_dict.values())
    output_dict[metadata_field] = list_of_kmer_values[int(kmer_row_index)]

#print json output
print(json.dumps(output_dict))

#save output to json file
with open(path_to_output, "wt") as f:
    json.dump(output_dict, f)