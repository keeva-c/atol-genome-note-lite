import json
from jinja2 import Template, Environment, FileSystemLoader, Undefined

print("Starting script")

# path_to_template = args.template
path_to_template = "dev/not-a-genome-note_template_scaffold-level.md"
path_to_genome_note_output = "results/genome_note_scaffold-level.md"
path_to_sample_metadata = "dev/sample_1.json"
path_to_supplement_template = "dev/hic_data_template.md"
path_to_supplement_output = "results/hic_data_for_genome_note.md"
path_to_supplement_metadata = "dev/sample_2.json"

# read the template
with open(path_to_template, "rt", encoding="utf-8") as f:
    template = Template(f.read())

# read the sample metadata
print(f"Reading sample metadata from {path_to_sample_metadata}")
with open(path_to_sample_metadata, "rt") as f:
    sample_metadata = json.load(f)

# render the template
print(f"Writing output to {path_to_genome_note_output}")
with open(path_to_genome_note_output, "wt", encoding="utf-8") as f:
    f.write(template.render(sample_metadata))

#repeat the above for suplementary metadata (for Hi-C data)

# read the template
with open(path_to_supplement_template, "rt", encoding="utf-8") as f:
    sup_template = Template(f.read())

# read the sample metadata
print(f"Reading supplementary sample metadata from {path_to_supplement_metadata}")
with open(path_to_supplement_metadata, "rt") as f:
    supplement_metadata = json.load(f)

# render the template
print(f"Writing supplementary output to {path_to_supplement_output}")
with open(path_to_supplement_output, "wt", encoding="utf-8") as f:
    f.write(sup_template.render(supplement_metadata))