# **A genome assembly for {{ "the " ~ organism.common_name ~ ", " if organism.common_name else "" }}*{{ organism.scientific_name }}*{{ " " ~ organism.authority if organism.authority else "" }}**

## **Authors**

{{ experiment.data_owner ~ ", " if experiment.data_owner else "" }}{{ sample.project_collaborators ~ ", " if sample.project_collaborators else "" }}Australian Tree of Life Infrastructure Capability, {{ sample.bpa_initiative }}

## **Abstract**

We have assembled and annotated a {% if assembly.scaffold_count==assembly.contig_count %}contig-level{% else %}scaffold-level{% endif %} genome
sequence for *{{ organism.scientific_name }}* ({{ organism.order_or_group ~ ": " if organism.order_or_group else "" }}{{
organism.family|default("",true)}}). The assembly comprises {% if assembly.scaffold_count!=assembly.contig_count %}{{ make_pretty_number(assembly.scaffold_count) }} scaffolds{% else %}{{ make_pretty_number(assembly.contig_count) }} contigs{% endif %} and spans {{ round_bases_up(assembly.genome_length) }}. It has a {% if assembly.scaffold_count!=assembly.contig_count %}scaffold N50 of {{ round_bases_up(assembly.scaffold_n50) }}, a {% endif %}contig N50 of {{ round_bases_up(assembly.contig_n50) }} and a BUSCO completeness score of {{ assembly.busco_c }}%. A total of {{ make_pretty_number(annotation.gene_count) }} genes were identified in annotation.

## **Introduction**

### **Species taxonomy**

{{ organism.tax_string|default("*higher taxon classification*",true) }}; *{{ organism.scientific_name }}*{{ " " ~ organism.authority if organism.authority else "" }} (NCBI:txid{{ organism.taxon_id|default("*ncbi taxon id*",true) }}).

### **Background**

The genome of {{ "the " ~ organism.common_name ~ ", " if organism.common_name else "" }}*{{ organism.scientific_name }}* was
sequenced as part of the {{ sample.bpa_initiative }} and has been assembled and annotated in collaboration with the Australian Tree of Life
Infrastructure Capability.

## **Genome sequence report**

Details about the assembled genome sequence, including key assembly metrics and target minimum standards set by the Earth Biogenome Project (EBP) (Earth BioGenome Project, 2024), are summarised in Table 1. {% if assembly.contact_map_image_path %}A Hi-C contact map for the assembly is provided in Figure 1.{% endif %} Genome annotation results are provided in Table 2.

| **Assembly information** | | **EBP standard** |
| --- | ---- | --- |
| Assembly name | {{ assembly.assembly_name }} | |
| Assembly accession | {{ assembly.assembly_accession }} | |
| Alternate haplotype assembly accession | {{ assembly.alt_hap_accession }} | |
| Span | {{ round_bases_up(assembly.genome_length) }} | |
| Number of gaps | {{ make_pretty_number(assembly.gap_count) }} | |
| Depth of coverage | {% if assembly.coverage %}{{ round_decimal(assembly.coverage) }}x{% else %}*not provided*{% endif %} | |
| Number of contigs | {{ make_pretty_number(assembly.contig_count) }} | |
| Contig N50 length | {{ round_bases_up(assembly.contig_n50) }} | > 1 Mb |
| Longest contig | {{ round_bases_up(assembly.longest_contig) }} | {% if assembly.scaffold_count!=assembly.contig_count %}
| Number of scaffolds | {{ make_pretty_number(assembly.scaffold_count) }} | |
| Scaffold N50 length | {{ round_bases_up(assembly.scaffold_n50) }} | Chromosomal scale |
| Longest scaffold | {{ round_bases_up(assembly.longest_scaffold) }} | {% endif %} | |
| Consensus quality (QV) | Primary assembly: {{ assembly.primary_qv }} | > 40 |
| | Alternate assembly: {{ assembly.alt_qv }} | |
| | Combined: {{ assembly.combined_qv }} | |
| *k*-mer completeness | Primary assembly: {{ assembly.primary_kmer }}% | > 90% |
| | Alternate assembly: {{ assembly.alt_kmer }}% | |
| | Combined: {{ assembly.combined_kmer }}% | |
| Organelles | {% if assembly.mito_size %}Mitochondrial genome: {{ make_pretty_number(assembly.mito_size) }} bp{% endif %}{% if assembly.plastid_size %} Plastid genome: {{ assembly.plastid_size }}{% endif %}{% if not assembly.mito_size and not assembly.plastid_size %}No organelles assembled{% endif %} | Complete single alleles |
| **BUSCO assessment** | |
| Full summary | {{ add_spaces(assembly.busco_string) }} |
| Single copy BUSCOs | {{ assembly.busco_single }}% | > 90% |
| Duplicated BUSCOs | {{ assembly.busco_duplicated }}% | < 5% |
| Reference set | {{ assembly.busco_ref }} | |
Table: Table 1: Genome assembly information for {{ assembly.assembly_name|default("*genome assembly name*",true) }}, sequenced from *{{ organism.scientific_name }}*.

{% if assembly.contact_map_image_path %}
![Figure 1: Hi-C contact map of the genome assembly, visualised using HiGlass.
Chromosomes are shown in order of size from left to right and top to bottom.]({{ assembly.contact_map_image_path }})
{% endif %}

| **Annotation information** | |
| -- | --- |
| Annotation name | {{ annotation.annotation_name }} |
| Number of genes | {{ make_pretty_number(annotation.gene_count) }} |
| Number of CDSs | {{ make_pretty_number(annotation.cds_count )}} |
| Number of gene transcripts | {{ make_pretty_number(annotation.transcript_count) }} |
| Average transcript length | {{ round_decimal(annotation.mean_transcript_length) }} |
| Average transcripts per gene | {{ round_decimal(annotation.mean_transcripts_per_gene) }} |
| Average exons per transcript | {{ round_decimal(annotation.mean_exons_per_transcript) }} |
| **BUSCO completeness** | |
| Full summary | {{ add_spaces(annotation.annot_busco_summary) }} |
| | *C: complete, S: single copy, D: duplicated, F: fragmented, M: missing* |
| Reference set | {{ annotation.annot_busco_lineage }} |
| **OMArk completeness** | |
| Full summary | {{ add_spaces(annotation.omark_completeness_summary) }} |
| | *S: single copy, D: duplicated [U: unexpected, E:expected], M: missing* |
| Number of HOGs | {{ make_pretty_number(annotation.conserved_hogs) }} |
| Reference lineage | {{ annotation.omark_lineage}} |
| **OMArk consistency** | |
| Number of proteins | {{ make_pretty_number(annotation.omark_protein_count) }}
| Consistent lineage placements | {{ round_decimal(annotation.omark_percent_consistent) }}% |
| Inconsistent lineage placements | {{ round_decimal(annotation.omark_percent_inconsistent) }}% |
| Contaminants | {{ round_decimal(annotation.omark_percent_contaminant) }}% |
| Unknown | {{ round_decimal(annotation.omark_percent_unknown) }}% |
Table: Table 2: Genome annotation information for {{ assembly.assembly_name|default("*genome assembly name*",true) }}.

## **Methods**

Information relating to sample collection, nucleic acid extraction, and
sequencing are provided in Tables 3, 4, and 5 respectively. {% if assembly.genomescope_image_path %}A frequency distribution graph of *k*-mers generated during sequencing is included in{% if assembly.contact_map_image_path %} Figure 2. {% else %} Figure 1. {% endif %}{% endif %}
An overview of the computational pipelines used in genome assembly, annotation and quality assessment are given in Table 6. Individual tools are listed in Table 7.

### **Sample acquisition**

| **Sample information** | |
| --- | ----- |
| **Sample: {{ sample.biosample_accession }}** | |
| BioSample accession | {{ sample.biosample_accession }} |
| Scientific name | *{{ organism.scientific_name }}* |
| TOLID | {{ sample.tolid }} |
| Specimen identifier | {{ sample.specimen_id }} |
| Specimen identifier defined by | {{ sample.specimen_id_description }} |
| Sex | {{ sample.sex }} |
| Lifestage | {{ sample.lifestage }} |
| Date | {{ sample.collection_date }}|
| Locality | {{ sample.region_and_locality }}{% if sample.region_and_locality!=sample.state_or_region %} |
| State/region | {{ sample.state_or_region }}{% else %}{% endif %} |
| Country/ocean | {{ sample.country_or_sea }} |
| Traditional place name | {{ sample.indigenous_location }} |
| Latitude | {{ sample.latitude }} |
| Longitude | {{ sample.longitude }} |
| Collected by | {{ sample.collected_by }} |
| Collection method | {{ sample.collection_method }} |
| Collection permit | {{ sample.collection_permit }} |
| Identified by | {{ sample.identified_by }} |
| Preservation method | {{ sample.preservation_method }} |
| Preservation temperature | {{ sample.preservation_temperature }} | 
{% include "supplementary_sample_data_for_genome_note.md" ignore missing %}
Table: Table 3: Sample information about the material used to generate
sequencing data.

### **Nucleic acid extraction**

| **Extraction information** | |
| --- | ----- |
| **Library: {{ experiment.bpa_library_id }}** | |
| Data type generated | {{ experiment.library_strategy }} |
| BioSample accession | {{ sample.biosample_accession }} |
| Sample tissue | {{ sample.organism_part}} |
| Nucleic acid extraction method | {{ sample.extraction_method }} |
| Nucleic acid treatment | {{ sample.nucleic_acid_treatment }} |
| Extract concentration (ng/ul) | {{ sample.nucleic_acid_conc }} |
{% include "supplementary_extract_data_for_genome_note.md" ignore missing %}
Table: Table 4: Methodological information about nucleic acid material
extracted for sequencing.

### **Sequencing**

| **Sequencing information** | |
| --- | ----- |
| **Run: {{ runs.sra_run_accession }}** | |
| Data type generated | {{ experiment.library_strategy }} |
| BioSample accession | {{ sample.biosample_accession }} |
| Library identifier | {{ experiment.bpa_library_id }} |
| Run accession | {{ runs.sra_run_accession }} |
| Read count | {% if runs.run_read_count %}{{ make_pretty_number(runs.run_read_count) }}{% else %}*not provided*{% endif %} |
| Base count | {% if runs.run_base_count %}{{ round_bases_up(runs.run_base_count) }}{% else %}*not provided*{% endif %} |
| Sequencing platform | {{ experiment.platform }} |
| Sequencing instrument | {{ experiment.instrument_model }} |
| Sequencing chemistry | {{ experiment.sequencing_kit }} |
| Sequencing facility | {{ experiment.GAL }} |
| Library preparation method | {{ experiment.library_construction_protocol }} |
| Library selection | {{ experiment.library_selection }} |
| Library layout | {{ experiment.library_layout }} |
| Library insert or fragment size | {{ experiment.insert_size }} |
| Flowcell type | {{ experiment.flowcell_type }} |
| Base caller model | {{ experiment.base_caller_model }} |
{% include "supplementary_seq_data_for_genome_note.md" ignore missing %}
Table: Table 5: Methodological information about sequencing runs.

{% if assembly.assembly.genomescope_image_path %}
{% if assembly.contact_map_image_path %}
![Figure 2: GenomeScope2.0 profile of *k*-mers generated during sequencing. Estimates of genome size (len), proportion of the genome sequence which is not repeated (uniq), and heterozygosity (ab), are given above the *k*-mer frequency distribution plot.]({{ assembly.genomescope_image_path }})
{% else %}
![Figure 1: GenomeScope2.0 profile of *k*-mers generated during sequencing. Estimates of genome size (len), proportion of the genome sequence which is not repeated (uniq), and heterozygosity (ab), are given above the *k*-mer frequency distribution plot.]({{ assembly.genomescope_image_path }})
{% endif %}
{% endif %}

### **Genome assembly and annotation**

|**Pipeline information** | |
| - | -- |
| **Genome assembly** | |
| - Pipeline | sanger-tol/genomeassembly |
| - Version | {{ assembly.assembly_pipeline_ver }} |
| - Source | [https://github.com/TomHarrop/atol_test_assembly](https://github.com/TomHarrop/atol_test_assembly) |
| - Adapted from | [https://github.com/sanger-tol/genomeassembly](https://github.com/sanger-tol/genomeassembly) |
| **Genome annotation** | |
| - Pipeline | {{ annotation.annotation_pipeline }} |
| - Version | {{ annotation.annotation_pipeline_ver }} |
| - Source | {{ annotation.annotation_pipeline_link }} |
| - Adapted from | tbd |
| **Annotation QC** | |
| - Pipeline | {{ annotation.annot_qc_pipeline }} |
| - Version | {{ annotation.annot_qc_pipeline_ver }} |
| - Source | {{ annotation.annot_qc_pipeline_link }} |
| - Adapted from | tbd |
Table: Table 6: Pipelines used in genome assembly, annotation and quality assessment.

| **Tool** | **Source** | **Reference** |
| -- | --- | --- |
| BEDTools | [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2) | Quinlan and Hall, 2010 |
| BUSCO | [https://gitlab.com/ezlab/busco](https://gitlab.com/ezlab/busco) | Manni *et al.*, 2021 |
| bwa-mem2 | [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | Vasimuddin *et al.*, 2019 |
| Cooler | [https://github.com/open2c/cooler](https://github.com/open2c/cooler) | Abdennur and Mirny, 2020 |
| FastK | [https://github.com/thegenemyers/FASTK](https://github.com/thegenemyers/FASTK) | NA |
| gawk | [https://www.gnu.org/software/gawk/](https://www.gnu.org/software/gawk/) | NA |
| GeneScopeFK | [https://github.com/thegenemyers/GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK) | NA |
| GenomeScope2.0 | [https://github.com/tbenavi1/genomescope2.0](https://github.com/tbenavi1/genomescope2.0) | Ranallo-Benavidez *et al.*, 2020 |
| Gfastats | [https://github.com/vgl-hub/gfastats](https://github.com/vgl-hub/gfastats) | Formenti *et al.*, 2022 |
| Hifiasm | [https://github.com/chhylp123/hifiasm](https://github.com/chhylp123/hifiasm) | Cheng *et al.*. 2021 |
| Juicer | [https://github.com/aidenlab/juicer](https://github.com/aidenlab/juicer) | Durand *et al.*, 2016 |
| JuicerTools | [https://github.com/aidenlab/JuicerTools](https://github.com/aidenlab/JuicerTools) | Durand *et al.*, 2016 |
| Merqury.FK | [https://github.com/thegenemyers/MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK) | Rhie *et al.*, 2020 |
| Minimap2 | [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2) | Li, 2018 |
| MitoFinder | [https://github.com/RemiAllio/MitoFinder](https://github.com/RemiAllio/MitoFinder) | Allio *et al.*, 2020 |
| MitoHiFi | [https://github.com/marcelauliano/MitoHiFi](https://github.com/marcelauliano/MitoHiFi) | Uliano-Silva *et al.*, 2023 |
| Nextflow | [https://github.com/nextflow-io/nextflow](https://github.com/nextflow-io/nextflow) | Di Tommaso *et al.*, 2017 |
| Nf-core | [https://nf-co.re/](https://nf-co.re/) | Ewels *et al.*, 2020 |
| Oatk | [https://github.com/c-zhou/oatk](https://github.com/c-zhou/oatk) | Zhou *et al.*, 2024 |
| pigz | [https://github.com/madler/pigz](https://github.com/madler/pigz) | NA |
| PretextMap | [https://github.com/sanger-tol/PretextMap](https://github.com/sanger-tol/PretextMap) | NA |
| PretextSnapshot | [https://github.com/sanger-tol/PretextSnapshot](https://github.com/sanger-tol/PretextSnapshot) | NA |
| purge_dups | [https://github.com/dfguan/purge_dups](https://github.com/dfguan/purge_dups) | Guan *et al.*, 2020 |
| samtools | [https://github.com/samtools/samtools](https://github.com/samtools/samtools) | Danecek *et al.*, 2021 |
| YaHS | [https://github.com/c-zhou/yahs](https://github.com/c-zhou/yahs) | Zhou *et al.*, 2023 |
Table: Table 9: Resources and software tools used in assembly pipelines.

## **Data availability**

Raw sequencing data, sample metadata, and genome assembly sequences are available from the European Nucleotide Archive under
the BioProject accession number {{ bioproject_accession }};
[https://identifiers.org/ena.embl/](https://identifiers.org/ena.embl/){{
bioproject_accession }}. Assembly and raw data accession identifiers are
reported in Tables 1 and 5. The genome annotation file is available from *TBD*. Raw sequence data and sample metadata were
originally submitted to the Bioplatforms Australia Data Portal
([https://data.bioplatforms.com/](https://data.bioplatforms.com/)),
and are available under the following data package identifiers: {{
experiment.bpa_package_id }}{% include "supplementary_package_data_for_genome_note.md" ignore missing %}.

The genome sequence is released openly for reuse.

## **Acknowledgements and funding information**

{% if sample.bpa_initiative=='Threatened Species Initiative' %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS) in partnership with the University of Sydney; the Australian Government Department of Climate Change, Energy, the Environment and Water; WA Department of Biodiversity, Conservation & Attractions; and Amazon Web Services.{% else %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS).{% endif %}

The genome has been assembled, annotated and published as part of the Australian Tree of Life Informatics Capacity, a platform provided by the Australian BioCommons.

Sample collection and preparation were supported by the individual project partners. The pipelines used to assemble and publish the genome sequence have been adapted from the original digital infrastructure developed as part of the Darwin Tree of Life project (Darwin Tree of Life Project Consortium, 2022). The generation of this sequence report has leveraged assets from the Tree of Life Genome Note pipeline (Babirye *et al.*, 2025).

## **References**

Abdennur, N. and Mirny, L. A. (2020) Cooler: Scalable storage for Hi-C
data and other genomically labeled arrays, *Bioinformatics*, 36 (1), pp.
311--316. DOI:10.1093/bioinformatics/btz540.

Allio, R., Schomaker‐Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz,
B. and Delsuc, F. (2020) MitoFinder: Efficient automated large‐scale
extraction of mitogenomic data in target enrichment phylogenomics,
*Molecular Ecology Resources*, 20 (4), pp. 892--905.
DOI:10.1111/1755-0998.13160.

Babirye, S. R., Chafin, T., Duong, C., Muffato, M., Qi, G., Sadasivan Baby, C., Surana, P., and Yates, B. (2025) sanger-tol/genomenote v2.1.1 (2.1.1). Zenodo. DOI:10.5281/zenodo.15052341.

Cheng, H., Concepcion, G. T., Feng, X., Zhang, H. and Li, H. (2021)
Haplotype-resolved *de novo* assembly using phased assembly graphs with
hifiasm, *Nature Methods*, 18 (2), pp. 170--175.
DOI:10.1038/s41592-020-01056-5.

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V.,
Pollard, M. O., *et al.* (2021) Twelve years of SAMtools and BCFtools,
*GigaScience*, 10 (2). DOI:10.1093/gigascience/giab008.

Darwin Tree of Life Project Consortium (2022) Sequence locally, think globally: The Darwin Tree of Life Project. *Proceedings of the National Academy of Sciences of the United States of America*, 119 (4), e2115642118. DOI:10.1073/pnas.2115642118.

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E.
and Notredame, C. (2017) Nextflow enables reproducible computational
workflows, *Nature Biotechnology*, 35 (4), pp. 316--319.
DOI:10.1038/nbt.3820.

Earth BioGenome Project (2024) Report on Assembly Standards Version 6.0 - September 2024. Retrieved 19 November, 2025, from [https://www.earthbiogenome.org/report-on-assembly-standards](https://www.earthbiogenome.org/report-on-assembly-standards)

Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm,
A., Garcia, M. U., Di Tommaso, P. and Nahnsen, S. (2020) The nf-core
framework for community-curated bioinformatics pipelines, *Nature
Biotechnology*, 38 (3), pp. 276--278. DOI:10.1038/s41587-020-0439-x.

Formenti, G., Abueg, L., Brajuka, A., Brajuka, N., Gallardo-Alba, C.,
Giani, A., Fedrigo, O. and Jarvis, E. D. (2022) Gfastats: conversion,
evaluation and manipulation of genome sequences using assembly graphs,
*Bioinformatics*, 38 (17), pp. 4214--4216.
DOI:10.1093/bioinformatics/btac460.

Guan, D., McCarthy, S. A., Wood, J., Howe, K., Wang, Y. and Durbin, R.
(2020) Identifying and removing haplotypic duplication in primary genome
assemblies., *Bioinformatics (Oxford, England)*, 36 (9), pp. 2896--2898.
DOI:10.1093/bioinformatics/btaa025.

Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences,
*Bioinformatics*, 34 (18), pp. 3094--3100.
DOI:10.1093/bioinformatics/bty191.

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A. and Zdobnov, E. M.
(2021) BUSCO update: Novel and streamlined workflows along with broader
and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic,
and viral genomes, *Molecular Biology and Evolution*, 38 (10), pp.
4647--4654. DOI:10.1093/molbev/msab199.

Quinlan, A. R. and Hall, I. M. (2010) BEDTools: a flexible suite of
utilities for comparing genomic features, *Bioinformatics*, 26 (6), pp.
841--842. DOI:10.1093/bioinformatics/btq033.

Ranallo-Benavidez, T. R., Jaron, K. S., & Schatz, M. C. (2020) GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature communications*, 11 (1), 1432. DOI:10.1038/s41467-020-14998-3.

Rhie, A., Walenz, B. P., Koren, S. and Phillippy, A. M. (2020) Merqury:
Reference-free quality, completeness, and phasing assessment for genome
assemblies, *Genome Biology*, 21 (1). DOI:10.1186/s13059-020-02134-9.

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