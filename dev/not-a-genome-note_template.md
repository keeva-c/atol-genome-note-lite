# **A genome assembly for the {{ organism.common_name }}, *{{ organism.scientific_name }}* {{ organism.authority }}**

## **Authors**

{{ project_lead }}, {{ project_collaborators }}, Australian Tree of Life
Infrastructure Capability, {{ bpa_initiative }}

## **Abstract**

We have assembled and annotated a {{ assembly_level }} level genome sequence
for the *{{ organism.scientific_name }}* ({{ organism.order_or_group }}: {{
organism.family }})*.* The assembly is comprised of ({{ contig_count }} contigs
or {{ scaffold_count }} scaffolds or {{ chromosome_count }} chromosomes) and
spans {{ genome_size }}Gb. It has a (scaffold N50 of {{ scaffold_n50 }} or
contig N50 of {{ contig_n50 }}) and a BUSCO completeness score of {{ busco_c
}}. We identified {{ pcg }} protein-coding genes and {{ ncg }} non-coding genes
in annotation.

# **Introduction**

## **Species taxonomy**

{{ tax_string }}; *{{ scientific_name }}* {{ authority }} (NCBI:txid{{
ncbi_taxid }}).

## **Background**

The genome of the {{ common_name }}, *{{ scientific_name }}*, was sequenced as
part of the {{ bpa_initiative }} and has been assembled and annotated in
collaboration with the Australian Tree of Life infrastructure capability.

# **Genome sequence report**

Details about the assembled genome sequence, including key assembly and
annotation metrics, are summarised in Table 1. A Hi-C contact map for the
assembly is provided in Figure 2.

Table 1: Genome assembly and annotation information for *{{ scientific_name
}}*, {{ assembly_name }}.

**Genome assembly** 

Assembly name\
 {{ assembly_name }}

Assembly accession\
 {{ assembly_accession }}

Alternate haplotype assembly accession\
 {{ alt_hap_accession }}

Span (Mb)\
 {{ genome_length }} 

Number of gaps\
 {{ gap_count }}

Depth of coverage\
 {{ coverage }}

Number of contigs\
 {{ contig_count }}

Contig N50 length (Mb)\
 {{ contig_n50 }}

Longest contig (Mb)\
 {{ longest_contig }}

Number of scaffolds\
 {{ scaffold_count }}

Scaffold N50 length (Mb)\
 {{ scaffold_n50 }}

Longest scaffold (Mb)\
 {{ longest_scaffold }}

Number of non-organelle chromosomes\
 {{ chromosome_count }}

**Assembly metrics\***

Consensus quality (QV)\
 {{ qv }}\
 *benchmark: 50*

*k*-mer completeness\
 {{ kmer }}%\
 *benchmark: 95%*

BUSCO\
 {{ busco_string }}\
 *benchmark: C = 95%*

BUSCO reference set\
 {{ busco_ref }}

Percentage of assembly mapped to chromosomes\
 {{ perc_assem }}%\
 *benchmark: 95%*

Sex chromosomes\
 {{ sex_chromosomes }}\
 *benchmark: localised homologous pairs*

Organelles\
 Mitochondrial genome: {{ mito_size }}kb\
 Plastid genome: {{ plastid_size }}kb\
 *or*\
 No organelles assembled\
 *benchmark: complete single alleles*

**Genome annotation**

Annotation name\
 {{ annotation_name }}

Number of protein-coding genes\
 {{ pcg }}

Number of non-coding genes\
 {{ ncg }}

Number of gene transcripts\
 {{ gene_transcripts }}
 
Average transcript length\
 {{ cds_length }}

Average number of coding transcripts per gene\
 {{ cds_per_gene }}

Average number of exons per transcript\
 {{ exon_per_transcript }}

\* Assembly metric benchmarks are adapted from column VGP-2020 of "Table 1:
Proposed standards and metrics for defining genome assembly quality" from Rhie
*et al.* (2021).

{{ contact_map_image }}

Figure 1: Hi-C contact map of the genome assembly, visualised using HiGlass.
Chromosomes are shown in order of size from left to right and top to bottom.

# **Methods**

Information relating to sample collection, nucleic acid extraction, and
sequencing are provided in Tables 2, 3, and 4 respectively. An overview of the
computational pipelines and workflows used in genome assembly and annotation
are given in Table 5. Individual tools are listed in Table 6.

## **Sample acquisition**

Table 2: Sample information about the material used to generate sequencing
data.

**Sample information**

Data type generated\
 WGS

Scientific name\
 *{{ scientific_name }}* 

TOLID\
 {{ tolid }}

BioSample accession\
 {{ biosample_accession }}

Specimen identifier\
 {{ specimen_voucher }}

Institution\
 {{ voucher_institution }}

Sex\
 {{ sex }}

Lifestage\
 {{ lifestage }}

Data type generated\
 Hi-C

Scientific name\
 *{{ scientific_name }}* 

TOLID\
 {{ tolid }}

BioSample accession\
 {{ biosample_accession }}

Specimen identifier\
 {{ specimen_voucher }}

Institution\
 {{ voucher_institution }}

Sex\
 {{ sex }}

Lifestage\
 {{ lifestage }}

Data type generated\
 RNA

Scientific name\
 *{{ scientific_name }}* 

TOLID\
 {{ tolid }}

BioSample accession\
 {{ biosample_accession }}

Specimen identifier\
 {{ specimen_voucher }}

Institution\
 {{ voucher_institution }}

Sex\
 {{ sex }}

Lifestage\
 {{ lifestage }}

**Collection information**  

Data type generated\
 WGS

Date\
 {{ collection_date }}
