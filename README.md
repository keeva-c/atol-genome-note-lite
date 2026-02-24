# AToL Genome Note Lite

Generates a markdown report for a genome sequence, including metadata and metrics about sampling, sequencing, assembly and annotation. This program is based on the [`sanger-tol/genomenote`](https://pipelines.tol.sanger.ac.uk/genomenote) pipeline developed by the Tree of Life informatics team at the Wellcome Sanger Institute.

## Usage

```
src/atol_genome_note_lite/explore.py 
   --wgs_metadata : 1 or more paths/to/wgs/metadata.json 
   --hic_metadata : 0 or more paths/to/hic/metadata.json 
   --output : optional path/to/results/output.md
   --w_annotation : specifies generating a genome note lite with annotation section 
   or 
   --wo_annotation : specifies generating a genome note lite without annotation section
```

To run [`explore.py`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/src/atol_genome_note_lite/explore.py) you need:
 * at least one JSON file containing metadata for the organism, sample, WGS run, assembly, and, optionally, annotation, formatted according to the AToL metadata schema

Optionally, you can include:
 * additional JSON files containing metadata for the organism, sample, and WGS run, if more than one set of reads was used to generate the assembly
 * one or more JSON files containing metadata for the organism, sample, and Hi-C run, if one or more sets of Hi-C reads were used to generate the assembly
 * specification of whether the genome note lite should include a section on genome annotation (if neither `--w_annotation` or `--wo_annotation` flags are used, no annotation information will be included by default)

## Output

 * a markdown report with key information about how the genome sequence was generated using descriptive text and tabular metadata and assembly metrics

 ## Full usage
 ```
usage: explore.py [-h] [--wgs_metadata WGS_METADATA [WGS_METADATA ...]] [--hic_metadata [HIC_METADATA ...]] [--output OUTPUT] [--w_annotation] [--wo_annotation]

This tool generates a draft genome note lite markdown document based on sample, read, assembly, and annotation metadata.

options:
  -h, --help            show this help message and exit
  --w_annotation        runs the genome note lite for an assembly with an annotation. (default: False)
  --wo_annotation       runs the genome note lite for an assembly only (no annotation). (default: False)

Input:
  --wgs_metadata WGS_METADATA [WGS_METADATA ...]
                        at least one JSON file listing metadata for a WGS sample and sequencing run. The first file should also contain assembly metadata and optionally annotation metadata. (default: None)
  --hic_metadata [HIC_METADATA ...]
                        optional JSON file/s listing metadata for a Hi-C sample and sequencing run. (default: None)

Output:
  --output OUTPUT       optional output file address. (default: results/genome_note_lite.md)
 ```

## Extracting file paths

The [`file_extractor.py`](x) script can be used to find the paths of the files used in the genome note lite given a text file listing the contents of the assembly pipeline output directory. Extracted paths only include outputs generated for hap1. Paths can be used in the `combined_parser.py` script to extract the relevant metrics or information.

> [!Note]
> This script can only retrieve file paths from output directories where the hi-c phasing option was selected during genome assembly.

### Full usage:

```
usage: file_extractor.py [-h] [--file_dir FILE_DIR] [--tolid TOLID] [--output OUTPUT]

This tool extracts file paths for files used in the genome note lite from a text file containing the directory contents generated in the genome assembly pipeline.

options:
  -h, --help           show this help message and exit
  --file_dir FILE_DIR  a text file listing all file paths contained in the genome assembly pipeline ouptut directory (default: None)
  --tolid TOLID        the ToLID for the specimen used to generate the genome assembly (default: None)
  --output OUTPUT      a JSON file containing the list of extracted paths with their file types (default: results/found_files.json)
```

## Assembly metadata and metrics

The [`combined_parser.py`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/src/atol_genome_note_lite/combined_parser.py) script can be used to parse ouput files genreated during the assembly and QC processes to extract and format the metrics and information required to generate the genome note lite. Assembly output files which can be parsed include:
 - BUSCO summaries, 
 - kmer stats,
 - QV stats, 
 - general assembly summaries, 
 - software version logs, 
 - contact maps, and
 - mitogenome stats.

### Full usage:

```
usage: combined_parser.py [-h] [--busco BUSCO] [--kmer KMER] [--qv QV] [--summary SUMMARY] [--software SOFTWARE] [--map [MAP]] [--kmer_plot [KMER_PLOT]] [--mito [MITO]] [--metadata [METADATA]] [--output OUTPUT] [--full_json FULL_JSON]

This tool parses files that are generated during the genome assembly pipeline and outputs them in a json format suitable for appending to an organism-sample-experiment-reads metadata json which can be used an input to the explore.py script

options:
  -h, --help            show this help message and exit

Input:
  --busco BUSCO         the JSON summary file generated during BUSCO analysis of the assembly - the file should end with busco.json (default: None)
  --kmer KMER           the kmer count file generated during merqury.FK analysis of the assembly - the file should end with .completeness.stats (default: None)
  --qv QV               the QV scores summary file generated during merqury.FK analysis of the assembly - the file should end with .qv (default: None)
  --summary SUMMARY     the summary text file generated when using the sanger_tol assembly pipeline - the file should end with .assembly_summary (default: None)
  --software SOFTWARE   the YAML file summarising the tools and versions used by the assembly pipeline - the file should end with genomeassembly_software_versions.yml (default: None)
  --map [MAP]           the PreText Snapshot PNG file of the contact map generated during assembly scaffolding - the file should end with FullMap.png (default: None)
  --kmer_plot [KMER_PLOT]
                        the GenomeScope2.0 PNG file containing a linear frequency distribution graph of kmers - the file should end with linear_plot.png (default: None)
  --mito [MITO]         the mitochondrion statistics TSV generated by MitoHiFi during mitogenome assembly - the file should end with contigs_stats.tsv (default: None)
  --metadata [METADATA]
                        the sample metadata JSON to which the assembly stats and metrics will be appended onto - should already include metadata for organism, sample, experiment, and runs (default: None)

Output:
  --output OUTPUT       a JSON output file including all assembly stats and information for inclusion in the genome note lite (default: results/assembly_metrics.json)
  --full_json FULL_JSON
                        the full JSON metadata object, including metadata for organism, sample, experiment, and runs, and assembly metrics (only generated if a metadata file is included as input) (default: results/full_metadata.json)
```

## How it works

The Genome Note Lite uses [Jinja](https://github.com/pallets/jinja) to populate a markdown template ([`templates/genome-note-lite-annot-template.md`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/templates/genome-note-lite-annot-template.md) or [`templates/genome-note-lite-asm-only-template.md`](https://github.com/keeva-c/atol-genome-note-lite/blob/main/templates/genome-note-lite-asm-only-template.md)) with variables specified in JSON input files. The template document is loosely based on the [`sanger-tol/genomenote`](https://pipelines.tol.sanger.ac.uk/genomenote) template, but has been reduced and re-formatted to better suit the requirements of the AToL project and its metadata. It is populated with the sample, sequencing, and assembly metadata, and calculated assembly metrics that are provided as input variables in the form of key-value pairs. Unlike the `sanger-tol/genomenote` pipeline, the Genome Note Lite does not run any bioinformatic analyses.

The input variables are preprocessed. This includes removing key-value pairs if the value is an empty string and adding an additional key-value pair with the full Bioplatforms Initiative name, if applicable. 

The template is populated with the metadata values provided for the primary WGS data package. If additional WGS data and/or Hi-C data are specified, the Genome Note Lite will populate supplementary templates using these metadata. The supplementary information is then included in the relevant sections of the main Genome Note Lite template. The supplementary markdown files are only generated and included where necessary. If any data were generated from the same sample, the sample metadata will not be replicated. If the additional data were generated from the same library, the nucleic acid extraction metadata will not be replicated.

During rendering, the template implements basic logic to determine certain wording choices. These largely depend on whether certain key-value pairs have been provided or whether the assembly is a contig- or scaffold-level assembly. Certain variables in the template have placeholder default values specified if the variable is not provided as input. If it is not specified in the template, the global default value for missing key-value pairs renders as "*not provided*". Additional formatting functions are also applied to standardise numerical representations.

After the template has been rendered, the Genome Note Lite deletes intermediate files including processed metadata and supplementary rendered templates.

The Genome Note Lite generates a single markdown file. If no file path is specified in the input arguments, the output will default to a file called `genome_note_lite.md` inside a `results/` directory. To convert the Genome Note to a `.docx` or `.pdf` format, you will need to use [Pandoc](https://pandoc.org/) or a similar conversion tool.