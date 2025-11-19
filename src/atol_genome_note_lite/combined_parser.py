#!/usr/bin/env python3

import json
import csv
import yaml

#set variables

path_to_busco_field_mapping = "dev/busco_to_fields.csv"
path_to_busco_file =  #input: busco.json file
busco_field_mapping_dict = {}
busco_results = {}
busco_output_dict = {}
path_to_busco_output = "dev/assembly_busco_fields.json"

path_to_kmer_field_mapping = "dev/kmer_to_fields.csv"
path_to_kmer_file =  #input: .completeness.stats file
kmer_field_mapping_dict = {}
kmer_values_dict = {}
kmer_output_dict = {}
path_to_kmer_output = "dev/assembly_kmer_fields.json"

path_to_qv_field_mapping = "dev/qv_to_fields.csv"
path_to_qv_file =  #input: .qv file
qv_field_mapping_dict = {}
qv_values_dict = {}
qv_output_dict = {}
path_to_qv_output = "dev/assembly_qv_fields.json"

path_to_summary_field_mapping = "dev/summary_to_fields.csv"
path_to_summary_file =  #input: .assembly_summary file
sep = ": "
summary_values_dict = {}
summary_output_dict = {}
summary_field_mapping_dict = {}
path_to_summary_output = "dev/assembly_summary_fields.json"

path_to_tools_file =  #input: genomeassembly_software_versions.yml file
workflow_version = {}
assembly_wf_ver_key = 'sanger-tol/genomeassembly'
path_to_tools_output = "dev/assembly_version_fields.json"

contact_map_file_path =  #input: FullMap.png file (if not applicable, set variable to None)

json_assembly_object = {}

#parse busco

#save genome note field names to a dictionary
with open(path_to_busco_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        busco_field_mapping_dict[line[1]]=line[0]

#extract "results" and "lineage_dataset" objects from json file and merge into one dictionary
with open(path_to_busco_file, "rt") as f:
    busco_stats = json.load(f)
    busco_results = busco_stats['results']
    busco_lineage = busco_stats['lineage_dataset']
    busco_for_mapping = busco_results | busco_lineage

#map values from merged dictionary to genome note field names
for metadata_field,busco_field_name in busco_field_mapping_dict.items():
    busco_output_dict[metadata_field] = busco_for_mapping[busco_field_name]

#add busco stats to json assembly object
json_assembly_object.update(busco_output_dict)

#parse kmer completeness

#save genome note field names to a dictionary
with open(path_to_kmer_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        kmer_field_mapping_dict[line[0]]=line[1]

#create a dictionary of kmer completeness values for primary, alt and combined assemblies
with open(path_to_kmer_file, "rt") as f:
    kmer_table = csv.reader(f, delimiter='\t')
    header_row = next(kmer_table)
    kmer_position = header_row.index("% Covered") #define the position of the kmer completeness values in each row
    for row in kmer_table:
        key = row[0]
        value = row[kmer_position]
        kmer_values_dict[key]=value

#map values from kmer completeness dictionary to genome note field names
for kmer_row_index,metadata_field in kmer_field_mapping_dict.items():
    list_of_kmer_values = list(kmer_values_dict.values())
    kmer_output_dict[metadata_field] = list_of_kmer_values[int(kmer_row_index)]

#add kmer values to json assembly object
json_assembly_object.update(kmer_output_dict)

#parse qv scores

#save genome note field names to a dictionary
with open(path_to_qv_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        qv_field_mapping_dict[line[1]]=line[0]

#create a dictionary of qv values for primary, alt and combined assemblies
with open(path_to_qv_file, "rt") as f:
    qv_table = csv.reader(f, delimiter='\t')
    header_row = next(qv_table)
    qv_position = header_row.index("QV") #define the position of the QV values in each row
    for row in qv_table:
        key = row[0]
        value = row[qv_position]
        qv_values_dict[key]=value

#map values from qv dictionary to genome note field names
for metadata_field,qv_row_index in qv_field_mapping_dict.items():
    list_of_qv_values = list(qv_values_dict.values())
    qv_output_dict[metadata_field] = list_of_qv_values[int(qv_row_index)]

#add qv values to json assembly object
json_assembly_object.update(qv_output_dict)

#parse summary

with open(path_to_summary_field_mapping, "rt") as f:
    csvreader = csv.reader(f)
    next(csvreader) #take out the header
    for line in csvreader:
        summary_field_mapping_dict[line[1]]=line[0]

with open(path_to_summary_file, "rt") as f:
    next(f) #take out the header
    for line in f:
        splits = line.strip().split(sep)
        if len(splits) == 2:
            key = splits[0]
            value = splits[1]
            summary_values_dict[key]=value

for metadata_field,summary_field_name in summary_field_mapping_dict.items():
    summary_output_dict[metadata_field] = summary_values_dict[summary_field_name]

#add assembly summary stats to json assembly object
json_assembly_object.update(summary_output_dict)

#parse tools

with open(path_to_tools_file, "rt") as f:
    tool_versions = yaml.load(f, Loader=yaml.SafeLoader)
    all_wf_versions = tool_versions['Workflow'] #extract all workflow-relevant versions in a dictionary

assembly_wf_ver = all_wf_versions[assembly_wf_ver_key] #extract value for assembly workflow version

workflow_version['assembly_pipeline_ver'] = assembly_wf_ver #create dictionary mapping genome note field name to assembly workflow version

#add workflow version to json assembly object
json_assembly_object.update(workflow_version)

#map contact map path
output_contact_map = {"contact_map_image_path": contact_map_file_path}

#add contact map path to json assembly object
json_assembly_object.update(output_contact_map)

#print json assembly object
print('"assembly":', json.dumps(json_assembly_object))