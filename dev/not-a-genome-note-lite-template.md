# **A genome assembly for {{ "the " ~ organism.common_name ~ ", " if organism.common_name else "" }}*{{ organism.scientific_name }}*{{ " " ~ organism.authority if organism.authority else "" }}**

## **Authors**

{{ experiment.data_owner ~ ", " if experiment.data_owner else "" }}{{ sample.project_collaborators ~ ", " if sample.project_collaborators else "" }}Australian Tree of Life
Infrastructure Capability, {{ sample.bpa_initiative }} Consortium

## **Abstract**

We have assembled a {% if assembly.scaffold_count==assembly.contig_count %}contig-level{% else %}scaffold-level{% endif %} genome
sequence for *{{ organism.scientific_name }}* ({{ organism.order_or_group ~ ": " if organism.order_or_group else "" }}{{
organism.family|default("",true)}}). The assembly comprises {% if assembly.scaffold_count!=assembly.contig_count %}{{ make_pretty_number(assembly.scaffold_count) }} scaffolds{% else %}{{ make_pretty_number(assembly.contig_count) }} contigs{% endif %} and spans {{ round_bases_up(assembly.genome_length) }}. It has a {% if assembly.scaffold_count!=assembly.contig_count %}scaffold N50 of {{ round_bases_up(assembly.scaffold_n50) }}, a {% endif %}contig N50 of {{ round_bases_up(assembly.contig_n50) }} and a BUSCO completeness score of {{ assembly.busco_c }}%.

# **Introduction**

## **Species taxonomy**

{{ organism.tax_string|default("*higher taxon classification*",true) }}; *{{ organism.scientific_name }}*{{ " " ~ organism.authority if organism.authority else "" }} (NCBI:txid{{
organism.taxon_id|default("*ncbi taxon id*",true) }}).

## **Background**

The genome of {{ "the " ~ organism.common_name ~ ", " if organism.common_name else "" }}*{{ organism.scientific_name }}*, was
sequenced as part of the {{ sample.bpa_initiative }} Consortium and has been assembled in collaboration with the Australian Tree of Life
Infrastructure Capability.

# **Genome sequence report**

Details about the assembled genome sequence, including key assembly 
metrics, are summarised in Table 1. {% if assembly.contact_map_image_path %}A Hi-C contact map for the assembly is provided in Figure 1.{% endif %}

Table 1: Genome assembly information for {{ assembly.assembly_name|default("*genome assembly name*",true) }}, sequenced from *{{ organism.scientific_name }}*.

**Genome assembly** 

Assembly name\
 {{ assembly.assembly_name }}

Assembly accession\
 {{ assembly.assembly_accession }}

Alternate haplotype assembly accession\
 {{ assembly.alt_hap_accession }}

Span\
 {{ round_bases_up(assembly.genome_length) }} 

Number of gaps\
 {{ make_pretty_number(assembly.gap_count) }}

Depth of coverage\
 {% if assembly.coverage %}{{ round_decimal(assembly.coverage) }}x{% else %}*not provided*{% endif %}

Number of contigs\
 {{ make_pretty_number(assembly.contig_count) }}

Contig N50 length\
 {{ round_bases_up(assembly.contig_n50) }}

Longest contig\
 {{ round_bases_up(assembly.longest_contig) }}
{% if assembly.scaffold_count!=assembly.contig_count %}
Number of scaffolds\
 {{ make_pretty_number(assembly.scaffold_count) }}

Scaffold N50 length\
 {{ round_bases_up(assembly.scaffold_n50) }}

Longest scaffold\
 {{ round_bases_up(assembly.longest_scaffold) }}
{% endif %}
**Assembly metrics\***

Consensus quality (QV)\
 Primary assembly: {{ assembly.primary_qv }}\
 Alternate assembly: {{ assembly.alt_qv }}\
 Combined: {{ assembly.combined_qv }}\
 *benchmark: 50*

*k*-mer completeness\
 Primary assembly: {{ assembly.primary_kmer }}%\
 Alternate assembly: {{ assembly.alt_kmer }}%\
 Combined: {{ assembly.combined_kmer }}%\
 *benchmark: 95%*

BUSCO\
 {{ assembly.busco_string }}\
 *benchmark: C = 95%*

BUSCO reference set\
 {{ assembly.busco_ref }}

Organelles\
 {% if assembly.mito_size %}Mitochondrial genome: {{ assembly.mito_size }}\{% endif %}{% if assembly.plastid_size %}Plastid genome: {{ assembly.plastid_size }}\{% endif %}{% if not assembly.mito_size and not assembly.plastid_size %}No organelles assembled\{% endif %}
 *benchmark: complete single alleles*

\* Assembly metric benchmarks are adapted from column VGP-2020 of "Table
1: Proposed standards and metrics for defining genome assembly quality"
from Rhie *et al.* (2021).
{% if assembly.contact_map_image_path %}
![image]({{ assembly.contact_map_image_path }})
Figure 1: Hi-C contact map of the genome assembly, visualised using HiGlass.
Chromosomes are shown in order of size from left to right and top to bottom.
{% endif %}

# **Methods**

Information relating to sample collection, nucleic acid extraction, and
sequencing are provided in Tables 2, 3, and 4 respectively. An overview
of the computational pipelines and workflows used in genome assembly are given in Table 5. Individual tools are listed in Table 6.

## **Sample acquisition**

Table 2: Sample information about the material used to generate
sequencing data.

**Sample information**

Data type generated\
 {{ experiment.library_strategy }}

Scientific name\
 *{{ organism.scientific_name }}* 

TOLID\
 {{ organism.tolid }}

BioSample accession\
 {{ sample.biosample_accession }}

Specimen identifier\
 {{ sample.specimen_voucher }}

Institution\
 {{ sample.voucher_institution }}

Sex\
 {{ sample.sex }}

Lifestage\
 {{ sample.lifestage }}

**Collection information**  

Date\
 {{ sample.collection_date }}

Location\
 Locality: {{ sample.region_and_locality }}\{% if sample.region_and_locality!=sample.state_or_region %}
 State/region: {{ sample.state_or_region }}\{% else %}{% endif %}
 Country/ocean: {{ sample.country_or_sea }}

Traditional place name\
 {{ sample.indigenous_location }}

Location coordinates\
 Latitude: {{ sample.latitude }}\
 Longitude: {{ sample.longitude }}

Collected by\
 {{ sample.collected_by }}

Collection method\
 {{ sample.collection_method }}

Collection permit\
 {{ sample.collection_permit }}

**Identification information**

Identified by\
 {{ sample.identified_by }}

**Preservation information**

Preservation method\
 {{ sample.preservation_method }}

Preservation temperature\
 {{ sample.preservation_temperature }}

{% include "supplementary_sample_data_for_genome_note.md" ignore missing %}

## **Nucleic acid extraction**

Table 3: Methodological information about nucleic acid material
extracted for sequencing.

Data type generated\
 {{ experiment.library_strategy }}

Sample tissue\
 {{ sample.organism_part}}

Nucleic acid extraction method\
 {{ sample.extraction_method }}

Nucleic acid treatment\
{{ sample.nucleic_acid_treatment }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

{% include "supplementary_extract_data_for_genome_note.md" ignore missing %}

## **Sequencing**

Table 4: Methodological information about sequencing runs.

Data type generated\
 {{ experiment.library_strategy }}

Run accession\
 {{ runs.sra_run_accession }}

Read count\
 {% if runs.run_read_count %}{{ make_pretty_number(runs.run_read_count) }}{% else %}*not provided*{% endif %}

Base count\
 {% if runs.run_base_count %}{{ round_bases_up(runs.run_base_count) }}{% else %}*not provided*{% endif %}

Sequencing instrument\
 {{ experiment.instrument_model }}

Sequencing chemistry\
 {{ experiment.sequencing_kit }}

Sequencing facility\
 {{ experiment.GAL }}

Library preparation method\
 {{ experiment.library_construction_protocol }}

Library selection\
 {{ experiment.library_selection }}

Library layout\
 {{ experiment.library_layout }}

Library insert or fragment size\
 {{ experiment.insert_size }}

Flowcell type\
 {{ experiment.flowcell_type }}

Base caller model\
 {{ experiment.base_caller_model }}

{% include "supplementary_seq_data_for_genome_note.md" ignore missing %}

## **Genome assembly**

Table 5: Pipelines and workflows used in genome assembly.

**Genome assembly**

Pipeline:
 sanger-tol/genomeassembly

Version:
 {{ assembly.assembly_pipeline_ver }}

Source:
 [https://github.com/TomHarrop/atol_test_assembly](https://github.com/TomHarrop/atol_test_assembly)

Adapted from:
 [https://github.com/sanger-tol/genomeassembly](https://github.com/sanger-tol/genomeassembly)


Table 6: Resources and software tools used in assembly pipelines.

Tool:
 BEDTools

Source: 
 [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)

Reference:
 Quinlan and Hall, 2010

Tool:
 BUSCO

Source:
 [https://gitlab.com/ezlab/busco](https://gitlab.com/ezlab/busco)

Reference:
 Manni *et al.*, 2021

Tool:
 bwa-mem2

Source:
 [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)

Reference:
 Vasimuddin *et al.*, 2019

Tool:
 Cooler

Source:
 [https://github.com/open2c/cooler](https://github.com/open2c/cooler)

Reference:
 Abdennur and Mirny, 2020

Tool:
 FastK

Source:
 [https://github.com/thegenemyers/FASTK](https://github.com/thegenemyers/FASTK)

Reference:
 NA

Tool:
 gawk

Source:
 [https://www.gnu.org/software/gawk/](https://www.gnu.org/software/gawk/)

Reference:
 NA

Tool:
 GeneScopeFK

Source:
 [https://github.com/thegenemyers/GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

Reference:
 NA

Tool:
 Gfastats

Source:
 [https://github.com/vgl-hub/gfastats](https://github.com/vgl-hub/gfastats)

Reference:
 Formenti *et al.*, 2022

Tool:
 Hifiasm

Source:
 [https://github.com/chhylp123/hifiasm](https://github.com/chhylp123/hifiasm)

Reference:
 Cheng *et al.*. 2021

Tool:
 Juicer

Source:
 [https://github.com/aidenlab/juicer](https://github.com/aidenlab/juicer)

Reference:
 Durand *et al.*, 2016

Tool:
 JuicerTools

Source:
 [https://github.com/aidenlab/JuicerTools](https://github.com/aidenlab/JuicerTools)

Reference:
 Durand *et al.*, 2016 

Tool:
 Merqury.FK

Source:
 [https://github.com/thegenemyers/MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

Reference:
 Rhie *et al.*, 2020

Tool:
 Minimap2

Source:
 [https://github.com/RemiAllio/MitoFinder](https://github.com/RemiAllio/MitoFinder)

Reference:
 Li, 2018

Tool:
 MitoFinder

Source:
 [https://github.com/RemiAllio/MitoFinder](https://github.com/RemiAllio/MitoFinder)

Reference:
 Allio *et al.*, 2020

Tool:
 MitoHiFi

Source:
 [https://github.com/marcelauliano/MitoHiFi](https://github.com/marcelauliano/MitoHiFi)

Reference:
 Uliano-Silva *et al.*, 2023

Tool:
 Nextflow

Source:
 [https://github.com/nextflow-io/nextflow](https://github.com/nextflow-io/nextflow)

Reference:
 Di Tommaso *et al.*, 2017

Tool:
 Nf-core

Source:
 [https://nf-co.re/](https://nf-co.re/)

Reference:
 Ewels *et al.*, 2020

Tool:
 Oatk

Source:
 [https://github.com/c-zhou/oatk](https://github.com/c-zhou/oatk)

Reference:
 Zhou *et al.*, 2024

Tool:
 pigz

Source:
 [https://github.com/madler/pigz](https://github.com/madler/pigz)

Reference:
 NA

Tool:
 PretextMap

Source:
 [https://github.com/sanger-tol/PretextMap](https://github.com/sanger-tol/PretextMap)

Reference:
 NA

Tool:
 PretextSnapshot

Source:
 [https://github.com/sanger-tol/PretextSnapshot](https://github.com/sanger-tol/PretextSnapshot)

Reference:
 NA

Tool:
 purge_dups

Source:
 [https://github.com/dfguan/purge_dups](https://github.com/dfguan/purge_dups)

Reference:
 Guan *et al.*, 2020

Tool:
 samtools

Source:
 [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

Reference:
 Danecek *et al.*, 2021

Tool:
 YaHS

Source:
 [https://github.com/c-zhou/yahs](https://github.com/c-zhou/yahs)

Reference:
 Zhou *et al.*, 2023

## **Data availability**

Raw sequencing data, sample metadata, and genome assembly sequences are available from the European Nucleotide Archive under
the BioProject accession number {{ bioproject_accession }};
[https://identifiers.org/ena.embl/](https://identifiers.org/ena.embl/){{
bioproject_accession }}. Assembly and raw data accession identifiers are
reported in Tables 1 and 4. Raw sequence data and sample metadata were
originally submitted to the Bioplatforms Australia Data Portal
([https://data.bioplatforms.com/](https://data.bioplatforms.com/)),
and are available under the following data package identifiers: {{
experiment.bpa_package_id }}{% include "supplementary_package_data_for_genome_note.md" ignore missing %}.

The genome sequence is released openly for reuse.

## **Acknowledgements and funding information**

{% if sample.bpa_initiative=='Threatened Species Initiative' %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} Consortium in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS) in partnership with the University of Sydney; the Australian Government Department of Climate Change, Energy, the Envrionment and Water; WA Department of Biodiversity, Conservation & Attractions; and Amazon Web Services.{% else %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} Consortium in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS).{% endif %}

The genome has been assembled and published as part of the Australian Tree of Life Informatics Capacity, a platform provided by the Australian BioCommons.

Sample collection and preparation were supported by the individual project partners. The pipelines used to assemble and publish the genome sequence have been adapted from the original digital infrastructure developed as part of the Darwin Tree of Life project (Darwin Tree of Life Project Consortium, 2022). The generation of this sequence report has leveraged assets from the Tree of Life Genome Note pipeline (Babirye *et al.*, 2025).

# **References**

Abdennur, N. and Mirny, L. A. (2020) Cooler: Scalable storage for Hi-C
data and other genomically labeled arrays, *Bioinformatics*, 36 (1), pp.
311--316. DOI:10.1093/bioinformatics/btz540.

Allio, R., Schomaker‐Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz,
B. and Delsuc, F. (2020) MitoFinder: Efficient automated large‐scale
extraction of mitogenomic data in target enrichment phylogenomics,
*Molecular Ecology Resources*, 20 (4), pp. 892--905.
DOI:10.1111/1755-0998.13160.

Altschul, S. F., Gish, W., Miller, W., Myers, E. W. and Lipman, D. J.
(1990) Basic local alignment search tool, *Journal of Molecular
Biology*, 215 (3), pp. 403--410. DOI:10.1016/S0022-2836(05)80360-2.

Babirye, S. R., Chafin, T., Duong, C., Muffato, M., Qi, G., Sadasivan Baby, C., Surana, P., and Yates, B. (2025) sanger-tol/genomenote v2.1.1 (2.1.1). Zenodo. DOI:10.5281/zenodo.15052341.

Bateman, A., Martin, M.-J., Orchard, S., Magrane, M., Ahmad, S., Alpi,
E., *et al.* (2023) UniProt: the Universal Protein Knowledgebase in
2023, *Nucleic Acids Research*, 51 (D1), pp. D523--D531.
DOI:10.1093/nar/gkac1052.

Buchfink, B., Reuter, K. and Drost, H.-G. (2021) Sensitive protein
alignments at tree-of-life scale using DIAMOND, *Nature Methods*, 18
(4), pp. 366--368. DOI:10.1038/s41592-021-01101-x.

Challis, R., Kumar, S., Sotero-Caio, C., Brown, M. and Blaxter, M.
(2023) Genomes on a Tree (GoaT): A versatile, scalable search engine for
genomic and sequencing project metadata across the eukaryotic tree of
life, *Wellcome Open Research*, 8, pp. 24.
DOI:10.12688/wellcomeopenres.18658.1.

Challis, R., Richards, E., Rajan, J., Cochrane, G. and Blaxter, M.
(2020) BlobToolKit -- interactive quality assessment of genome
assemblies, *G3: Genes, Genomes, Genetics*, 10 (4), pp. 1361--1374.
DOI:10.1534/g3.119.400908.

Cheng, H., Concepcion, G. T., Feng, X., Zhang, H. and Li, H. (2021)
Haplotype-resolved *de novo* assembly using phased assembly graphs with
hifiasm, *Nature Methods*, 18 (2), pp. 170--175.
DOI:10.1038/s41592-020-01056-5.

da Veiga Leprevost, F., Grüning, B. A., Alves Aflitos, S., Röst, H. L.,
Uszkoreit, J., Barsnes, H., *et al.* (2017) BioContainers: an
open-source and community-driven framework for software standardization,
*Bioinformatics*, 33 (16), pp. 2580--2582.
DOI:10.1093/bioinformatics/btx192.

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V.,
Pollard, M. O., *et al.* (2021) Twelve years of SAMtools and BCFtools,
*GigaScience*, 10 (2). DOI:10.1093/gigascience/giab008.

Darwin Tree of Life Project Consortium (2022) Sequence locally, think globally: The Darwin Tree of Life Project. *Proceedings of the National Academy of Sciences of the United States of America*, 119 (4), e2115642118. DOI:10.1073/pnas.2115642118.

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E.
and Notredame, C. (2017) Nextflow enables reproducible computational
workflows, *Nature Biotechnology*, 35 (4), pp. 316--319.
DOI:10.1038/nbt.3820.

Diesh, C., Stevens, G. J., Xie, P., De Jesus Martinez, T., Hershberg, E.
A., Leung, A., *et al.* (2023) JBrowse 2: a modular genome browser with
views of synteny and structural variation, *Genome Biology*, 24 (1), pp.
74. DOI:10.1186/s13059-023-02914-z.

Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm,
A., Garcia, M. U., Di Tommaso, P. and Nahnsen, S. (2020) The nf-core
framework for community-curated bioinformatics pipelines, *Nature
Biotechnology*, 38 (3), pp. 276--278. DOI:10.1038/s41587-020-0439-x.

Ewels, P., Magnusson, M., Lundin, S. and Käller, M. (2016) MultiQC:
summarize analysis results for multiple tools and samples in a single
report, *Bioinformatics*, 32 (19), pp. 3047--3048.
DOI:10.1093/bioinformatics/btw354.

Formenti, G., Abueg, L., Brajuka, A., Brajuka, N., Gallardo-Alba, C.,
Giani, A., Fedrigo, O. and Jarvis, E. D. (2022) Gfastats: conversion,
evaluation and manipulation of genome sequences using assembly graphs,
*Bioinformatics*, 38 (17), pp. 4214--4216.
DOI:10.1093/bioinformatics/btac460.

Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J.,
Tomkins-Tinch, C. H., Valieris, R. and Köster, J. (2018) Bioconda:
sustainable and comprehensive software distribution for the life
sciences, *Nature Methods*, 15 (7), pp. 475--476.
DOI:10.1038/s41592-018-0046-7.

Guan, D., McCarthy, S. A., Wood, J., Howe, K., Wang, Y. and Durbin, R.
(2020) Identifying and removing haplotypic duplication in primary genome
assemblies., *Bioinformatics (Oxford, England)*, 36 (9), pp. 2896--2898.
DOI:10.1093/bioinformatics/btaa025.

Harry, E. (2022) PretextView (Paired REad TEXTure Viewer): A desktop
application for viewing pretext contact maps. Available from:
[https://github.com/wtsi-hpag/PretextView](https://github.com/wtsi-hpag/PretextView).

Kerpedjiev, P., Abdennur, N., Lekschas, F., McCallum, C., Dinkla, K.,
Strobelt, H., *et al.* (2018) HiGlass: web-based visual exploration and
analysis of genome interaction maps, *Genome Biology*, 19 (1), pp. 125.
DOI:10.1186/s13059-018-1486-1.

Kurtzer, G. M., Sochat, V. and Bauer, M. W. (2017) Singularity:
Scientific containers for mobility of compute, *PLOS ONE*, 12 (5), pp.
e0177459. DOI:10.1371/journal.pone.0177459.

Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences,
*Bioinformatics*, 34 (18), pp. 3094--3100.
DOI:10.1093/bioinformatics/bty191.

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A. and Zdobnov, E. M.
(2021) BUSCO update: Novel and streamlined workflows along with broader
and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic,
and viral genomes, *Molecular Biology and Evolution*, 38 (10), pp.
4647--4654. DOI:10.1093/molbev/msab199.

Merkel, D. (2014) Docker: lightweight Linux containers for consistent
development and deployment, *Linux Journal*, 2014 (239), pp. 2. DOI:10.5555/2600239.2600241.

Pointon, D.-L., Eagles, W., Sims, Y., Muffato, M. and Surana, P. (2023)
sanger-tol/treeval v1.0.0 -- Ancient Atlantis.
DOI:10.5281/zenodo.10047653.

Quinlan, A. R. and Hall, I. M. (2010) BEDTools: a flexible suite of
utilities for comparing genomic features, *Bioinformatics*, 26 (6), pp.
841--842. DOI:10.1093/bioinformatics/btq033.

Ranallo-Benavidez, T. R., Jaron, K. S., & Schatz, M. C. (2020) GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature communications*, 11 (1), 1432. DOI:10.1038/s41467-020-14998-3.

Rhie, A., McCarthy, S. A., Fedrigo, O., Damas, J., Formenti, G., Koren,
S., *et al.* (2021) Towards complete and error-free genome assemblies of
all vertebrate species, *Nature*, 592 (7856), pp. 737--746.
DOI:10.1038/s41586-021-03451-0.

Rhie, A., Walenz, B. P., Koren, S. and Phillippy, A. M. (2020) Merqury:
Reference-free quality, completeness, and phasing assessment for genome
assemblies, *Genome Biology*, 21 (1). DOI:10.1186/s13059-020-02134-9.

Sayers, E. W., Cavanaugh, M., Clark, K., Pruitt, K. D., Sherry, S. T.,
Yankie, L. and Karsch-Mizrachi, I. (2024) GenBank 2024 Update, *Nucleic
Acids Research*, 52 (D1), pp. D134--D137. DOI:10.1093/nar/gkad903.

Sims, Y., Butt, Z., Chafin, T., Challis, R., Kumar, S., Muffato, M., Qi, G., Ramos Díaz, A., Surana, P., and Yates, B. (2025) sanger-tol/blobtoolkit v0.8.0 – Sprigatito (0.8.0). Zenodo. DOI:10.5281/zenodo.15466127.

Surana, P., Muffato, M. and Qi, G. (2023b) sanger-tol/readmapping:
sanger-tol/readmapping v1.1.0 - Hebridean Black (1.1.0). Zenodo.
DOI:10.5281/zenodo.7755669.

Uliano-Silva, M., Ferreira, J. G. R. N., Krasheninnikova, K., Blaxter,
M., Mieszkowska, N., Hall, N., *et al.* (2023) MitoHiFi: a python
pipeline for mitochondrial genome assembly from PacBio high fidelity
reads, *BMC Bioinformatics*, 24 (1), pp. 288.
DOI:10.1186/s12859-023-05385-y.

Vasimuddin, Md., Misra, S., Li, H. and Aluru, S. (2019) Efficient
Architecture-Aware Acceleration of BWA-MEM for Multicore Systems, in:
*2019 IEEE International Parallel and Distributed Processing Symposium
(IPDPS)*. IEEE, pp. 314--324.

Zhou, C., Brown, M., Blaxter, M., The Darwin Tree of Life Consortium, McCarthy, S. A. and Durbin, R. (2024) Oatk: a de novo assembly tool for complex plant organelle genomes, *bioRxiv*. DOI:10.1101/2024.10.23.619857

Zhou, C., McCarthy, S. A. and Durbin, R. (2023) YaHS: yet another Hi-C
scaffolding tool, *Bioinformatics*, 39 (1).
DOI:10.1093/bioinformatics/btac808.
