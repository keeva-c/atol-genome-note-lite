# **A genome assembly for {{ "the " ~ taxonomy_info.ncbi_common_name ~ ", " if taxonomy_info.ncbi_common_name else "" }}*{{ taxonomy_info.ncbi_scientific_name }}*{{ " " ~ taxonomy_info.ncbi_authority if taxonomy_info.ncbi_authority else "" }}**

## **Authors**

{{ experiment.data_owner ~ ", " if experiment.data_owner else "" }}{{ sample.project_collaborators ~ ", " if sample.project_collaborators else "" }}Australian Tree of Life Bioinformatics group, {{ sample.bpa_initiative }}

## **Abstract**

We have assembled a {{ assembly.assembly_level }}-level genome sequence for *{{ taxonomy_info.ncbi_scientific_name }}* ({{ taxonomy_info.ncbi_order }}:{{ taxonomy_info.ncbi_family }}). The assembly comprises {% if assembly.assembly_level=='scaffold' %}{{ make_pretty_number(assembly.scaffold_count) }} scaffolds{% elif assembly.assembly_level=='contig' %}{{ make_pretty_number(assembly.contig_count) }} contigs{% endif %} and spans {{ round_bases_up(assembly.genome_length) }}. It has a {% if assembly.assembly_level=='scaffold' %}scaffold N50 of {{ round_bases_up(assembly.scaffold_n50) }}, a {% endif %}contig N50 of {{ round_bases_up(assembly.contig_n50) }} and a BUSCO completeness score of {{ assembly.busco_c }}%.

## **Introduction**

### **Species taxonomy**

{{ taxonomy_info.ncbi_full_lineage }}; *{{ taxonomy_info.ncbi_scientific_name }}*{{ " " ~ taxonomy_info.ncbi_authority if taxonomy_info.ncbi_authority else "" }} (NCBI:txid{{ taxonomy_info.ncbi_taxon_id }}).

### **Background**

The genome of {{ "the " ~ taxonomy_info.ncbi_common_name ~ ", " if taxonomy_info.ncbi_common_name else "" }}*{{ taxonomy_info.ncbi_scientific_name }}* was
sequenced as part of the {{ sample.bpa_initiative }} and has been assembled in collaboration with the Australian Tree of Life
Bioinformatics group.{% if assembly.rna_data_available %} As part of this project, additional RNA-seq data were generated and made available to support downstream applications of the genome assembly. Details of how these data were generated are included in the Methods section.{% endif %}

## **Methods**

Information relating to sample collection, nucleic acid extraction and library preparation, and sequencing are provided in Tables 1, 2, and 3 respectively.{% if assembly.genomescope_image_path %} A frequency distribution graph of *k*-mers generated during sequencing is included in Figure 1.{% endif %} An overview of the computational pipelines used in genome assembly and quality assessment are given in Table 4.

### **Sample acquisition**

| **Sample information** | |
| --- | ----- |
| **Sample: {{ sample.biosample_accession }}** | |
| Scientific name | *{{ taxonomy_info.ncbi_scientific_name }}* |
| ToLID | {{ sample.tolid }} |
| Specimen identifier | {{ sample.specimen_id }} |
| Specimen identifier defined by | {{ sample.specimen_id_description }} |
| Specimen-level BioSample accession | {{ specimen.biosample_accession }} |
| Sample-level BioSample accession | {{ sample.biosample_accession }} |
| Sex | {{ sample.sex }} |
| Lifestage | {{ sample.lifestage }} |
| Collection date | {{ sample.collection_date }} |
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
| Sample tissue | {{ sample.organism_part}} |
{% include "supplementary_sample_data_for_genome_note.md" ignore missing %}
Table: Table 1: Sample information about the material used to generate
sequencing data.

### **Nucleic acid extraction and library preparation**

| **Extraction and library information** | |
| --- | ----- |
| **Library: {{ experiment.bpa_library_id }}** | |
| Data type generated | {{ experiment.library_strategy }} |
| Sample-level BioSample accession | {{ sample.biosample_accession }} |
| Nucleic acid extraction method | {{ sample.extraction_method }} |
| Nucleic acid treatment | {{ sample.nucleic_acid_treatment }} |
| Extract concentration (ng/ul) | {{ sample.nucleic_acid_conc }} |
| Library preparation method | {{ experiment.library_construction_protocol }} |
| Library selection method | {{ experiment.library_selection }} |
| Library insert or fragment size | {{ experiment.insert_size }} |
{% include "supplementary_extract_data_for_genome_note.md" ignore missing %}
Table: Table 2: Methodological information about nucleic acid material
extracted for sequencing and library preparation.

### **Sequencing**

| **Sequencing information** | |
| --- | ----- |
| **Run: {{ runs.sra_run_accession }}** | |
| Data type generated | {{ experiment.library_strategy }} |
| Sample-level BioSample accession | {{ sample.biosample_accession }} |
| Library identifier | {{ experiment.bpa_library_id }} |
| Run accession | {{ runs.sra_run_accession }} |
| Read count | {% if runs.run_read_count %}{{ make_pretty_number(runs.run_read_count) }}{% else %}*not provided*{% endif %} |
| Base count | {% if runs.run_base_count %}{{ round_bases_up(runs.run_base_count) }}{% else %}*not provided*{% endif %} |
| Sequencing platform | {{ experiment.platform }} |
| Sequencing instrument | {{ experiment.instrument_model }} |
| Sequencing chemistry | {{ experiment.sequencing_kit }} |
| Sequencing facility | {{ experiment.GAL }} |
| Library layout | {{ experiment.library_layout }} |
| Flowcell type | {{ experiment.flowcell_type }} | {% if experiment.platform=='Oxford Nanopore' %}
| Base caller model | {{ experiment.base_caller_model }} | {% endif %}
{% include "supplementary_seq_data_for_genome_note.md" ignore missing %}
Table: Table 3: Methodological information about sequencing runs.

{% if assembly.genomescope_image_path %}
![Figure 1: GenomeScope2.0 profile of *k*-mers generated during sequencing. Estimates of genome size (len), proportion of the genome sequence which is not repeated (uniq), and heterozygosity (ab), are given above the *k*-mer frequency distribution plot.]({{ assembly.genomescope_image_path }})
{% endif %}

### **Genome assembly**

|**Pipeline information** | |
| - | -- |
| **Reads QC** | | {% if experiment.platform=='PacBio' %}
| - Pipeline | amytims/atol-qc-raw-pacbio |
| - Version | *tbd* |
| - Source | [https://github.com/amytims/atol-qc-raw-pacbio](https://github.com/amytims/atol-qc-raw-pacbio) | {% elif experiment.platform=='Oxford Nanopore' %}
| - Pipeline | TomHarrop/atol-qc-raw-ont |
| - Version | *tbd* |
| - Source | [https://github.com/TomHarrop/atol-qc-raw-ont](https://github.com/TomHarrop/atol-qc-raw-ont) | {% endif %}{% if assembly.assembly_level!='contig' %}
| - Pipeline | TomHarrop/atol-qc-raw-shortread |
| - Version | *tbd* |
| - Source | [https://github.com/TomHarrop/atol-qc-raw-shortread](https://github.com/TomHarrop/atol-qc-raw-shortread) | {% endif %}
| **Genome assembly** | |
| - Pipeline | sanger-tol/genomeassembly |
| - Version | {{ assembly.genomeassembly_pipeline_version }} |
| - Reference | Downie *et al.*, 2025 |
| - Source | [https://github.com/sanger-tol/genomeassembly](https://github.com/sanger-tol/genomeassembly) |
| **Decontamination** | |
| - Pipeline | sanger-tol/ascc |
| - Version | {{ assembly.ascc_pipeline_version }} |
| - Reference | Pointon *et al.*, 2026 |
| - Source | [https://github.com/sanger-tol/ascc](https://github.com/sanger-tol/ascc) | {% if assembly.assembly_level!='contig' %}
| **Assembly visualisation** | |
| - Pipeline | sanger-tol/treeval |
| - Version | {{ assembly.treeval_pipeline_version }} |
| - Reference | Pointon *et al.*, 2026 |
| - Source | [https://github.com/sanger-tol/treeval](https://github.com/sanger-tol/treeval) | {% endif %}
Table: Table 4: Pipelines used in genome assembly and quality assessment.

## **Genome sequence report**

Details about the assembled genome sequence, including key assembly metrics and target minimum standards set by the Earth BioGenome Project (EBP) (Earth BioGenome Project, 2026), are summarised in Table 5. {% if assembly.contact_map_image_path %}A Hi-C contact map for the assembly is provided in {% if assembly.genomescope_image_path %}Figure 2.{% else %} Figure 1.{% endif %}{% endif %}

| **Assembly information** | | **EBP standard** |
| --- | ---- | -- |
| Assembly name | {{ assembly.assembly_name }} | |
| Assembly accession | {{ assembly.assembly_accession }} | |
| Alternate haplotype assembly accession | {{ assembly.alt_hap_accession }} | |
| Span | {{ round_bases_up(assembly.genome_length) }} | | {% if assembly.assembly_level=='scaffold' %}
| Number of gaps | {{ make_pretty_number(assembly.gap_count) }} | | {% endif %} 
| Depth of coverage | {% if assembly.coverage %}{{ round_decimal(assembly.coverage) }}x{% else %}*not provided*{% endif %} | |
| Number of contigs | {{ make_pretty_number(assembly.contig_count) }} | |
| Contig N50 length | {{ round_bases_up(assembly.contig_n50) }} | > 1 Mb |
| Longest contig | {{ round_bases_up(assembly.longest_contig) }} | | {% if assembly.assembly_level=='scaffold' %}
| Number of scaffolds | {{ make_pretty_number(assembly.scaffold_count) }} | |
| Scaffold N50 length | {{ round_bases_up(assembly.scaffold_n50) }} | Chromosomal scale |
| Longest scaffold | {{ round_bases_up(assembly.longest_scaffold) }} | | {% endif %} |
| Consensus quality (QV) | Primary assembly: {{ assembly.primary_qv }} | > 40 |
| | Alternate assembly: {{ assembly.alt_qv }} | |
| | Combined: {{ assembly.combined_qv }} | |
| *k*-mer completeness | Primary assembly: {{ assembly.primary_kmer }}% | > 90% |
| | Alternate assembly: {{ assembly.alt_kmer }}% | |
| | Combined: {{ assembly.combined_kmer }}% | |
| Organelles | {% if assembly.mito_size %}Mitochondrial genome: {{ make_pretty_number(assembly.mito_size) }} bp{% endif %}{% if assembly.plastid_size %} Plastid genome: {{ assembly.plastid_size }}{% endif %}{% if not assembly.mito_size and not assembly.plastid_size %}No organelles assembled{% endif %} | Complete single alleles |
| Full BUSCO summary | {{ add_spaces(assembly.busco_string) }} |
| | *C: complete, S: single copy, D: duplicated/multi-copy, F: fragmented, M: missing, n: number of markers, E: proportion with internal stop codons* | |
| Single copy BUSCOs | {{ assembly.busco_single }}% | > 90% |
| Duplicated BUSCOs | {{ assembly.busco_duplicated }}% | < 5% |
| BUSCO Reference set | {{ assembly.busco_ref }} | |
Table: Table 5: Genome assembly information for {{ assembly.assembly_name|default("*genome assembly name*",true) }}, sequenced from *{{ taxonomy_info.ncbi_scientific_name }}*.

{% if assembly.contact_map_image_path %}
{% if assembly.genomescope_image_path %}
![Figure 2: Hi-C contact map of the genome assembly, visualised using HiGlass.
{% if assembly.assembly_level=='chromosome' %}Chromosomes {% else %}Scaffolds {% endif %}are shown in order of size from left to right and top to bottom.]({{ assembly.contact_map_image_path }})
{% else %}
![Figure 1: Hi-C contact map of the genome assembly, visualised using HiGlass.
{% if assembly.assembly_level=='chromosome' %}Chromosomes {% else %}Scaffolds {% endif %}are shown in order of size from left to right and top to bottom.]({{ assembly.contact_map_image_path }})
{% endif %}
{% endif %}

## **Data and code availability**

Raw sequencing data, sample metadata, and genome assembly sequences are available from the European Nucleotide Archive under
the BioProject accession number {{ project.project_accession }};
[https://identifiers.org/ena.embl/](https://identifiers.org/ena.embl/){{
project.project_accession }}. Assembly and raw data accession identifiers are
reported in Tables 1, 3, and 5. Raw sequence data and sample metadata were
originally submitted to the Bioplatforms Australia Data Portal
([https://data.bioplatforms.com/](https://data.bioplatforms.com/)),
and are available under the following data package identifiers: {{
experiment.bpa_package_id }}{% include "supplementary_package_data_for_genome_note.md" ignore missing %}.

The genome sequence is released openly for reuse.

A code repository containing the workflow configuration files used to generate the assembly is hosted on GitHub at [https://github.com/AToL-Bioinformatics/{{ sample.tolid }}.{{ assembly.assembly_version }}](https://github.com/AToL-Bioinformatics/{{ sample.tolid }}.{{ assembly.assembly_version }}). Original bioinformatics piplines are available from the links provided in Table 4.

## **Acknowledgements and funding information**

{% if sample.bpa_initiative=='Threatened Species Initiative' %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS) in partnership with the University of Sydney; the Australian Government Department of Climate Change, Energy, the Environment and Water; WA Department of Biodiversity, Conservation & Attractions; and Amazon Web Services.{% else %}We would like to acknowledge the contribution of the {{ sample.bpa_initiative }} in the generation of data used in this publication. The Initiative is supported by funding from Bioplatforms Australia, enabled by the Commonwealth Government National Collaborative Research Infrastructure Strategy (NCRIS).{% endif %}

The genome has been assembled and published using infrastructure provided by the Australian Tree of Life Bioinformatics group, an initiative of the Australian BioCommons.

Sample collection and preparation were supported by the individual project partners. The pipelines used to assemble and publish the genome sequence have been adapted from the original digital infrastructure developed as part of the Darwin Tree of Life project (Darwin Tree of Life Project Consortium, 2022). The generation of this sequence report has leveraged assets from the Tree of Life Genome Note pipeline (Babirye *et al.*, 2025).

## **References**

Babirye, S. R., Chafin, T., Duong, C., Muffato, M., Qi, G., Sadasivan Baby, C., Surana, P., and Yates, B. (2025) sanger-tol/genomenote v2.1.1 (2.1.1). Zenodo. DOI:10.5281/zenodo.15052341.

Darwin Tree of Life Project Consortium (2022) Sequence locally, think globally: The Darwin Tree of Life Project. *Proceedings of the National Academy of Sciences of the United States of America*, 119 (4), e2115642118. DOI:10.1073/pnas.2115642118.

Downie, J., Krashennikova, K., Muffato, M., Qi, G., Sims, Y., & Surana, P. (2025). sanger-tol/genomeassembly v0.50.0 - Threadtail (0.50.0). Zenodo. DOI:10.5281/zenodo.10391851.

Earth BioGenome Project (2024) Report on Assembly Standards Version 6.0 - September 2024. Retrieved 19 November, 2025, from [https://www.earthbiogenome.org/report-on-assembly-standards](https://www.earthbiogenome.org/report-on-assembly-standards){% if assembly.assembly_level!='contig' %}

Damon-Lee Pointon, Will Eagles, YSims, Matthieu Muffato, Mahesh Binzer-Panchal, Priyanka Surana, Guoying Qi, & Tree of Life service account. (2026). sanger-tol/treeval: 1.4.7 - Ancient Hippaforalkus (H7) (1.4.7). DOI:10.5281/zenodo.10047653.{% endif %}

Damon-Lee Pointon, YSims, eeaunin, Jim Downie, & Will Eagles. (2026). sanger-tol/ascc: 1.6.0 - Red Microphone (0.6.0). Zenodo. DOI:10.5281/zenodo.14883765.