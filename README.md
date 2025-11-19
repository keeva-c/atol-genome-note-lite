# AToL Genome Note Lite

This repository contains the code to generate a simplified Genome Note.

The purpose of a Genome Note is to provide information about the sampling, sequencing, and assembly of a genome sequence in the form of a short publication. The Genome Note Lite program is based 
on the [sanger-tol/genomenote](https://pipelines.tol.sanger.ac.uk/genomenote) pipeline developed by the Tree of Life informatics team at the Wellcome Sanger Institute. The Genome Note Lite differs from the sanger-tol Genome Note pipeline in that it does not include any bioinformatics analyses on the assembly sequence or input data, instead populating a basic template document with key metadata and assembly statistics in a largely tabular format.

The Genome Note Lite is being developed as part of the [Australian Tree of Life](https://www.biocommons.org.au/atol) (AToL) - an infrastructure project aiming to expedite genome assembly, annotation, and publishing by bridging and automating processes in bioinformatics, data brokering, and publication.

This code remains under development.

## Overview

The Genome Note Lite uses [Jinja](https://github.com/pallets/jinja) to populate a markdown template ([`not-a-genome-note-lite-template.md`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/dev/not-a-genome-note-lite-template.md)) with variables specified in JSON input files. The template document is loosely based on the sanger-tol/genomenote template, but has been reduced and re-formatted to better suit the requirements of the AToL project and its metadata. It is populated with the sample, sequencing, and assembly metadata, and calculated assembly metrics that are provided as input variables in the form of key-value pairs. 

The input variables are preprocessed. This includes removing key-value pairs if the value is an empty string and adding an additional key-value pair with the full Bioplatforms Initiative name, if applicable. 

The template is populated with the metadata values provided for the primary WGS data package. If additional WGS data and/or Hi-C data are specified, the Genome Note Lite will populate supplementary templates based on these metadata. The supplementary information is then included in the relevant sections of the main Genome Note Lite template. The supplementary markdown files are only generated and included where necessary - if the additional WGS or Hi-C data were generated from the same sample as the primary WGS data, the sample metadata will not be replicated, if the additional data were generated from the same library, the nucleic acid extraction metadata will not be replicated.

During rendering, the template implements basic logic to determine certain wording choices. These largely depend on whether certain key-value pairs have been provided or whether the assembly is a contig- or scaffold-level assembly. Certain variables in the template have placeholder default values specified if the variable is not provided as input. If it is not specified in the template, the global default value for missing key-value pairs renders as "*not provided*". Additional formatting functions are also applied to standardise numerical representations.

After the template has been rendered, the Genome Note Lite wipes the supplementary markdown files to eliminate the risk that these are unintentionally included in subsequent runs.

## Input

The Genome Note Lite takes nested JSON files as input, including JSON objects for *Sample*, *Organism*, *Experiment*, *Reads*, *Assembly*, and *Annotation*. Each object contains key-value pairs for metadata and metrics relating to each category, as specified by the AToL Schema.

The paths to the JSON files are specified in the main script [`explore.py`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/src/atol_genome_note_lite/explore.py) as:
 - `path_to_sample_metadata` (mandatory) (the JSON containing information about the primary WGS sample used and the sequencing data generated from it)
 - `path_to_WGS_supplement_metadata` (optional) (the JSON containing information about the sample and sequencing data generated in addition to the primary WGS data)
 - `path_to_hic_supplement_metadata` (optional) (the JSON containing information about the sample and sequencing data generated for Hi-C)

If there is no supplemental WGS or Hi-C data, these variable must be set to `None`.

Assembly and annotation information are only taken from the file specified in `path_to_sample_metadata`.

### Assembly metadata and metrics

The parse scripts under `src/atol_genome_note_lite/` can be used to extract and format the relevant *Assembly* metadata and metrics from the output files generated from running the [sanger-tol/genomeassembly](https://pipelines.tol.sanger.ac.uk/genomeassembly) pipeline. The parse scripts have been combined in [`combined_parser.py`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/src/atol_genome_note_lite/combined_parser.py). File paths to the relevant assembly outputs need to be set as variables inside the script. The expected file extensions are specified in the comments. The CSV assets needed to map the data from the assembly pipeline output files to the AToL schema JSON are in the `dev/` directory.

## Output

The Genome Note Lite generates a single markdown file, `genome_note_lite.md` inside a `results/` directory. To convert the Genome Note to a `.docx` or `.pdf` format, you will need to use [Pandoc](https://pandoc.org/), or a similar conversion tool.