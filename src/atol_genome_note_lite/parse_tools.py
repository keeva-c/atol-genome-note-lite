#!/usr/bin/env python3

import yaml
import json

path_to_file = "dev/software_versions.yml"
workflow_version = {}
assembly_wf_ver_key = 'sanger-tol/genomeassembly'
path_to_output = "dev/assembly_version_fields.json"

with open(path_to_file, "rt") as f:
    tool_versions = yaml.load(f, Loader=yaml.SafeLoader)
    all_wf_versions = tool_versions['Workflow'] #extract all workflow-relevant versions in a dictionary

assembly_wf_ver = all_wf_versions[assembly_wf_ver_key] #extract value for assembly workflow version

workflow_version['assembly_pipeline_ver'] = assembly_wf_ver #create dictionary mapping genome note field name to assembly workflow version

#save output to json file
with open(path_to_output, "wt") as f:
    json.dump(workflow_version, f)