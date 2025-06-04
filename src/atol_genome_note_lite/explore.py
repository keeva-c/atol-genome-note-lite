import json
from jinja2 import Template

print("Starting script")

# path_to_template = args.template
path_to_template = "dev/not-a-genome-note_template.md"
path_to_output = "results/genome_note.md"
path_to_sample_metadata = "dev/sample_1.json"

# read the template
with open(path_to_template, "rt") as f:
    template = Template(f.read())

# read the sample metadata
print(f"Reading sample metadata from {path_to_sample_metadata}")
with open(path_to_sample_metadata, "rt") as f:
    sample_metadata = json.load(f)

# render the template
print(f"Writing output to {path_to_output}")
with open(path_to_output, "wt") as f:
    f.write(template.render(sample_metadata))
