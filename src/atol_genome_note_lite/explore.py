import json
from jinja2 import Template, Environment, FileSystemLoader, Undefined

print("Starting script")

#setting the global default undefined value
class undefined_tokens(Undefined):
    def __str__(self):
        return "*not provided*"

#functions to make things in the template pretty
def make_pretty_number(ugly_number):
    if float(ugly_number) > 1_000:
        return "{:,}".format(int(ugly_number))

def round_decimal(decimal):
    return "{:.1f}".format(float(decimal))

def round_bases_up(base_pairs):
    if int(base_pairs) > 1_000_000_000:
        return "{:.2f} Gb".format(int(base_pairs) / 1_000_000_000)
    elif int(base_pairs) > 1_000_000:
        return "{:.2f} Mb".format(int(base_pairs) / 1_000_000)
    elif int(base_pairs) > 1_000:
        return "{:.2f} kb".format(int(base_pairs) / 1_000)
    else:
        return base_pairs + " bp"

# path_to_template = args.template
path_to_supplement_metadata = "dev/sample_2.json"
path_to_sample_metadata = "dev/sample_1.json"
path_to_sample_supplement_output = "dev/hic_sample_data_for_genome_note.md"
path_to_extract_supplement_output = "dev/hic_extract_data_for_genome_note.md"
path_to_seq_supplement_output = "dev/hic_seq_data_for_genome_note.md"
path_to_genome_note_output = "results/genome_note_lite.md"

# setting the environment for the genome note templates
env = Environment(loader=FileSystemLoader("dev"),undefined=undefined_tokens)

# read the supplementary templates
sup_samp_template = env.get_template("sup_sample_template.md")
sup_ext_template = env.get_template("sup_extract_template.md")
sup_seq_template = env.get_template("sup_sequencing_template.md")

# read the supplementary sample metadata
print(f"Reading supplementary sample metadata from {path_to_supplement_metadata}")
with open(path_to_supplement_metadata, "rt") as f:
    supplement_metadata = json.load(f)

# render the supplementary templates
sup_sample_render = sup_samp_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
sup_ext_render = sup_ext_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)
sup_seq_render = sup_seq_template.render(supplement_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal)

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
    f.write(template.render(sample_metadata,make_pretty_number=make_pretty_number,round_bases_up=round_bases_up,round_decimal=round_decimal))