import json
from jinja2 import Template, Environment, FileSystemLoader

print("Starting script")

# path_to_template = args.template
path_to_sample_supplement_output = "dev/hic_sample_data_for_genome_note.md"
path_to_extract_supplement_output = "dev/hic_extract_data_for_genome_note.md"
path_to_seq_supplement_output = "dev/hic_seq_data_for_genome_note.md"
path_to_supplement_metadata = "dev/sample_2.json"
path_to_sample_metadata = "dev/sample_1.json"
path_to_genome_note_output = "results/genome_note_scaffold-level.md"

env = Environment(loader=FileSystemLoader("dev"))

# read the supplementary templates
sup_samp_template = env.get_template("sup_sample_template.md")
sup_ext_template = env.get_template("sup_extract_template.md")
sup_seq_template = env.get_template("sup_sequencing_template.md")

# read the supplementary sample metadata
print(f"Reading supplementary sample metadata from {path_to_supplement_metadata}")
with open(path_to_supplement_metadata, "rt") as f:
    supplement_metadata = json.load(f)

# render the supplementary templates
sup_sample_render = sup_samp_template.render(supplement_metadata)
sup_ext_render = sup_ext_template.render(supplement_metadata)
sup_seq_render = sup_seq_template.render(supplement_metadata)

print(f"Writing supplementary output to {path_to_sample_supplement_output}")
with open(path_to_sample_supplement_output, "wt", encoding="utf-8") as f:
    f.write(sup_sample_render)

print(f"Writing supplementary output to {path_to_extract_supplement_output}")
with open(path_to_extract_supplement_output, "wt", encoding="utf-8") as f:
    f.write(sup_ext_render)

print(f"Writing supplementary output to {path_to_seq_supplement_output}")
with open(path_to_seq_supplement_output, "wt", encoding="utf-8") as f:
    f.write(sup_seq_render)

# read the core template
template = env.get_template("not-a-genome-note_template_scaffold-level.md")

# read the sample metadata
print(f"Reading sample metadata from {path_to_sample_metadata}")
with open(path_to_sample_metadata, "rt") as f:
    sample_metadata = json.load(f)

# render the core template
print(f"Writing output to {path_to_genome_note_output}")
with open(path_to_genome_note_output, "wt", encoding="utf-8") as f:
    f.write(template.render(sample_metadata))