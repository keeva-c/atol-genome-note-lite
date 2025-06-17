# **A genome assembly for the {{ organism.common_name }}, *{{ organism.scientific_name }}* {{ organism.authority }}**

## **Authors**

{{ sample.project_lead }}, {{ sample.project_collaborators }}, Australian Tree of Life
Infrastructure Capability, {{ sample.bpa_initiative }}

## **Abstract**

We have assembled and annotated a {{ assembly.assembly_level }} level genome
sequence for the *{{ organism.scientific_name }}* ({{ organism.order_or_group }}: {{
organism.family }}). The assembly is comprised of ({{ assembly.contig_count }} contigs
or {{ assembly.scaffold_count }} scaffolds or {{ assembly.chromosome_count }} chromosomes)
and spans {{ assembly.genome_size }}Gb. It has a (scaffold N50 of {{ assembly.scaffold_n50
}} or contig N50 of {{ assembly.contig_n50 }}) and a BUSCO completeness score of
{{ assembly.busco_c }}. We identified {{ annotation.pcg }} protein-coding genes and {{ annotation.ncg
}} non-coding genes in annotation.

# **Introduction**

## **Species taxonomy**

{{ organism.tax_string }}; *{{ organism.scientific_name }}* {{ organism.authority }} (NCBI:txid{{
organism.ncbi_taxid }}).

## **Background**

The genome of the {{ organism.common_name }}, *{{ organism.scientific_name }}*, was
sequenced as part of the {{ sample.bpa_initiative }} project and has been assembled and
annotated in collaboration with the Australian Tree of Life
infrastructure capability.

# **Genome sequence report**

Details about the assembled genome sequence, including key assembly and
annotation metrics, are summarised in Table 1. A Hi-C contact map for
the assembly is provided in Figure 1.

Table 1: Genome assembly and annotation information for {{ assembly.assembly_name }}, sequenced from *{{ organism.scientific_name }}*.

**Genome assembly** 

Assembly name\
 {{ assembly.assembly_name }}

Assembly accession\
 {{ assembly.assembly_accession }}

Alternate haplotype assembly accession\
 {{ assembly.alt_hap_accession }}

Span (Mb)\
 {{ assembly.genome_length }} 

Number of gaps\
 {{ assembly.gap_count }}

Depth of coverage\
 {{ assembly.coverage }}

Number of contigs\
 {{ assembly.contig_count }}

Contig N50 length (Mb)\
 {{ assembly.contig_n50 }}

Longest contig (Mb)\
 {{ assembly.longest_contig }}

Number of scaffolds\
 {{ assembly.scaffold_count }}

Scaffold N50 length (Mb)\
 {{ assembly.scaffold_n50 }}

Longest scaffold (Mb)\
 {{ assembly.longest_scaffold }}

Number of non-organelle chromosomes\
 {{ assembly.chromosome_count }}

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

Percentage of assembly mapped to chromosomes\
 {{ assembly.perc_assem }}%\
 *benchmark: 95%*

Sex chromosomes\
 {{ assembly.sex_chromosomes }}\
 *benchmark: localised homologous pairs*

Organelles\
 Mitochondrial genome: {{ assembly.mito_size }}kb\
 Plastid genome: {{ assembly.plastid_size }}kb\
 *or*\
 No organelles assembled\
 *benchmark: complete single alleles*

**Genome annotation**

Annotation name\
 {{ annotation.annotation_name }}

Number of protein-coding genes\
 {{ annotation.pcg }}

Number of non-coding genes\
 {{ annotation.ncg }}

Number of gene transcripts\
 {{ annotation.gene_transcripts }}
 
Average transcript length\
 {{ annotation.cds_length }}

Average number of coding transcripts per gene\
 {{ annotation.cds_per_gene }}

Average number of exons per transcript\
 {{ annotation.exon_per_transcript }}

\* Assembly metric benchmarks are adapted from column VGP-2020 of "Table
1: Proposed standards and metrics for defining genome assembly quality"
from Rhie *et al.* (2021).

![image]({{ assembly.contact_map_image }})

Figure 1: Hi-C contact map of the genome assembly, visualised using HiGlass.
Chromosomes are shown in order of size from left to right and top to bottom.

# **Methods**

Information relating to sample collection, nucleic acid extraction, and
sequencing are provided in Tables 2, 3, and 4 respectively. An overview
of the computational pipelines and workflows used in genome assembly and
annotation are given in Table 5. Individual tools are listed in Table 6.

## **Sample acquisition**

Table 2: Sample information about the material used to generate
sequencing data.

**Sample information**

Data type generated\
 WGS

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

Data type generated\
 Hi-C

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

Data type generated\
 RNA

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

Data type generated\
 WGS

Date\
 {{ sample.collection_date }}

Location\
 {{ sample.region_and_locality }}, {{ sample.country_or_sea }}

Traditional place name\
 {{ sample.indigenous_location }}

Location coordinates\
 {{ sample.latitude }}, {{ sample.longitude }}

Collected by\
 {{ sample.collected_by }}

Collection method\
 {{ sample.collection_method }}

Collection permit\
 {{ sample.collection_permit }}

Data type generated\
 Hi-C

Date\
 {{ sample.collection_date }}

Location\
 {{ sample.region_and_locality }}, {{ sample.country_or_sea }}

Traditional place name\
 {{ sample.indigenous_location }}

Location coordinates\
 {{ sample.latitude }}, {{ sample.longitude }}

Collected by\
 {{ sample.collected_by }}

Collection method\
 {{ sample.collection_method }}

Collection permit\
 {{ sample.collection_permit }}

Data type generated\
 RNA

Date\
 {{ sample.collection_date }}

Location\
 {{ sample.region_and_locality }}, {{ sample.country_or_sea }}

Traditional place name\
 {{ sample.indigenous_location }}

Location coordinates\
 {{ sample.latitude }}, {{ sample.longitude }}

Collected by\
 {{ sample.collected_by }}

Collection method\
 {{ sample.collection_method }}

Collection permit\
 {{ sample.collection_permit }}

**Identification information**

Data type generated\
 WGS

Identified by\
 {{ sample.identified_by }}

Data type generated\
 Hi-C

Identified by\
 {{ sample.identified_by }}

Data type generated\
 RNA

Identified by\
 {{ sample.identified_by }}

**Preservation information**

Data type generated\
 WGS

Preservation method\
 {{ sample.preservation_method }}

Preservation temperature\
 {{ sample.preservation_temperature }}

Data type generated\
 Hi-C

Preservation method\
 {{ sample.preservation_method }}

Preservation temperature\
 {{ sample.preservation_temperature }}

Data type generated\
 RNA

Preservation method\
 {{ sample.preservation_method }}

Preservation temperature\
 {{ sample.preservation_temperature }}

## **Nucleic acid extraction**

*For older data retrofitted to the AToL schema:*

Table 3: Methodological information about nucleic acid material
extracted for sequencing.

Data type generated\
 WGS

Sample tissue\
 {{ sample.organism_part}}

Nucleic acid extraction method\
 {{ sample.extraction_method }}

Nucleic acid treatment\
{{ sample.nucleic_acid_treatment }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

Data type generated\
 Hi-C

Sample tissue\
 {{ sample.organism_part}}

Nucleic acid extraction method\
 {{ sample.extraction_method }}

Nucleic acid treatment\
{{ sample.nucleic_acid_treatment }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

Data type generated\
 RNA

Sample tissue\
 {{ sample.organism_part}}

Nucleic acid extraction method\
 {{ sample.extraction_method }}

Nucleic acid treatment\
{{ sample.nucleic_acid_treatment }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

*For data produced under the new BPA/AToL schema:*

Table 3: Methodological information about nucleic acid material
extracted for sequencing.

Data type generated\
 WGS

Sample tissue\
 {{ sample.organism_part }}

Nucleic acid extraction protocol\
 DOI:{{ sample.extraction_protocol_DOI }}

Extract volume (ul)\
 {{ sample.nucleic_acid_volume }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

Data type generated\
 Hi-C

Sample tissue\
 {{ sample.organism_part }}

Nucleic acid extraction protocol\
 DOI:{{ sample.extraction_protocol_DOI }}

Extract volume (ul)\
 {{ sample.nucleic_acid_volume }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

Data type generated\
 RNA

Sample tissue\
 {{ sample.organism_part }}

Nucleic acid extraction protocol\
 DOI:{{ sample.extraction_protocol_DOI }}

Extract volume (ul)\
 {{ sample.nucleic_acid_volume }}

Extract concentration (ng/ul)\
 {{ sample.nucleic_acid_conc }}

## **Sequencing**

Table 4: Methodological information about sequencing runs.

Data type generated\
 WGS

Run accession\
 {{ runs.sra_run_accession }}

Read count\
 {{ runs.run_read_count }}

Base count (Gb)\
 {{ runs.run_base_count}}

Sequencing instrument\
 {{ experiment.instrument_model }}

Sequencing chemistry\
 {{ experiment.sequencing_kit }}

Sequencing facility\
 {{ sample.GAL }}

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

Data type generated\
 Hi-C

Run accession\
 {{ runs.sra_run_accession }}

Read count\
 {{ runs.run_read_count }}

Base count (Gb)\
 {{ runs.run_base_count}}

Sequencing instrument\
 {{ experiment.instrument_model }}

Sequencing chemistry\
 {{ experiment.sequencing_kit }}

Sequencing facility\
 {{ sample.GAL }}

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

Data type generated\
 RNA

Run accession\
 {{ runs.sra_run_accession }}

Read count\
 {{ runs.run_read_count }}

Base count (Gb)\
 {{ runs.run_base_count}}

Sequencing instrument\
 {{ experiment.instrument_model }}

Sequencing chemistry\
 {{ experiment.sequencing_kit }}

Sequencing facility\
 {{ sample.GAL }}

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

## **Genome assembly and annotation**

Table 5: Pipelines and workflows used in genome assembly and annotation.

**Genome assembly**

Pipeline:
 {{ assembly.assembly_pipeline }}

Version:
 {{ assembly.assembly_pipeline_ver }}

Source:
 {{ assembly.assembly_pipeline_link }}

Adapted from:
 [https://github.com/sanger-tol/genomeassembly](https://github.com/sanger-tol/genomeassembly)

**Assembly decontamination**

Pipeline:
 {{ assembly.decontamination_pipeline }}

Version:
 {{ assembly.decontamination_pipeline_ver }}
 
Source:
{{ assembly.decontamination_pipeline_link }}

Adapted from:
 [https://github.com/sanger-tol/ascc](https://github.com/sanger-tol/ascc)

**Manual curation**

Workflow:
 {{ assembly.contact_map_pipeline }}

Version:
 {{ assembly.contact_map_pipeline_ver }}

Source:
 {{ assembly.contact_map_pipeline_link }}

Adapted from:
 [https://github.com/sanger-tol/PretextView](https://github.com/sanger-tol/PretextView) (Harry, 2022) 

Workflow:
 {{ assembly.curation_pipeline }}

Version:
 {{ assembly.curation_pipeline_ver }}

Source:
 {{ assembly.curation_pipeline_link }}

Adapted from:
 [https://gitlab.com/wtsi-grit/rapid-curation](https://gitlab.com/wtsi-grit/rapid-curation)

**Assembly evaluation**

Pipeline:
 {{ assembly.blob_tk_pipeline }}

Version:
 {{ assembly.blob_tk_pipeline_ver }}

Source:
 {{ assembly.blob_tk_pipeline_link }}

Adapted from:
 [https://github.com/blobtoolkit/blobtoolkit](https://github.com/blobtoolkit/blobtoolkit) (Challis *et al.*, 2020) and [https://pipelines.tol.sanger.ac.uk/blobtoolkit](https://pipelines.tol.sanger.ac.uk/blobtoolkit) (Sims *et al.*, 2025)

Pipeline:
 {{ assembly.read_mapping_pipeline }}

Version:
 {{ assembly.read_mapping_pipeline_ver }}

Source:
 {{ assembly.read_mapping_pipeline_link }}

Adapted from:
 [https://github.com/sanger-tol/readmapping](https://github.com/sanger-tol/readmapping) (Surana *et al.*, 2023b)

Pipeline:
 {{ assembly.treeval_pipeline }}

Version:
 {{ assembly.treeval_pipeline_ver }}

Source:
 {{ assembly.treeval_pipeline_link }}

Adapted from:
[https://github.com/sanger-tol/treeval](https://github.com/sanger-tol/treeval) (Pointon *et al.*, 2023)

**Genome annotation**

Workflow:
 {{ annotation.annotation_pipeline }}

Version:
 {{ annotation.annotation_pipeline_ver }}

Source:
 {{ annotation.annotation_pipeline_link }}

Adapted from:
 *TBD*

Table 6: Resources and software tools used in assembly and annotation
pipelines.

Tool:
 BEDTools

Version:
 {{ assembly.BEDTools_version}}

Source: 
 [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)

Reference:
 Quinlan and Hall, 2010

Tool:
 BioConda

Version:
 NA

Source: 
 [https://github.com/bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)

Reference:
 Grüning *et al.*, 2018

Tool:
 BioContainers

Version:
 NA

Source: 
 [https://biocontainers.pro/](https://biocontainers.pro/)

Reference:
 da Veiga Leprevost *et al.*, 2017

Tool:
 BLAST

Version:
 {{ assembly.BLAST_version}}

Source: 
 [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

Reference:
 Altschul *et al.*, 1990

Tool:
 BUSCO

Version:
 {{ assembly.BUSCO_version }}

Source:
 [https://gitlab.com/ezlab/busco](https://gitlab.com/ezlab/busco)

Reference:
 Manni *et al.*, 2021

Tool:
 bwa-mem2

Version:
 {{ assembly.bwamem2_version }}

Source:
 [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)

Reference:
 Vasimuddin *et al.*, 2019

Tool:
 Cooler

Version:
 {{ assembly.Cooler_version }}

Source:
 [https://github.com/open2c/cooler](https://github.com/open2c/cooler)

Reference:
 Abdennur and Mirny, 2020

Tool:
 DIAMOND

Version:
 {{ assembly.DIAMOND_version }}

Source:
 [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)

Reference:
 Buchfink *et al.*, 2021

Tool:
Docker

Version:
 {{ assembly.Docker_version }}

Source:
 [https://www.docker.com/](https://www.docker.com/)

Reference:
 Merkel, 2014

Tool:
fasta_windows

Version:
 {{ assembly.fasta_windows_version }}

Source:
 [https://github.com/tolkit/fasta_windows](https://github.com/tolkit/fasta_windows)

Reference:
 NA

Tool:
 FastK

Version:
 {{ assembly.FastK_version }}

Source:
 [https://github.com/thegenemyers/FASTK](https://github.com/thegenemyers/FASTK)

Reference:
 NA

Tool:
 gawk

Version:
 {{ assembly.gawk_version }}

Source:
 [https://www.gnu.org/software/gawk/](https://www.gnu.org/software/gawk/)

Reference:
 NA

Tool:
 GeneScopeFK

Version:
 {{ assembly.GeneScopeFK_version }}

Source:
 [https://github.com/thegenemyers/GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)

Reference:
 NA

Tool:
 GenomeScope2

Version:
 {{ assembly.GenomeScope2_version }}

Source:
 [https://github.com/tbenavi1/genomescope2.0](https://github.com/tbenavi1/genomescope2.0)

Reference:
 Ranallo-Benavidez *et al.*, 2020

Tool:
 Gfastats

Version:
 {{ assembly.Gfastats_version }}

Source:
 [https://github.com/vgl-hub/gfastats](https://github.com/vgl-hub/gfastats)

Reference:
 Formenti *et al.*, 2022

Tool:
 GoaT CLI

Version:
 {{ assembly.GoaT_CLI_version }}

Source:
 [https://github.com/genomehubs/goat-cli](https://github.com/genomehubs/goat-cli)

Reference:
 Challis *et al.*, 2023

Tool:
 Hifiasm

Version:
 {{ assembly.Hifiasm_version }}

Source:
 [https://github.com/chhylp123/hifiasm](https://github.com/chhylp123/hifiasm)

Reference:
 Cheng *et al.*. 2021

Tool:
 HiGlass

Version:
 {{ assembly.HiGlass_version }}

Source:
 [https://github.com/higlass/higlass](https://github.com/higlass/higlass)

Reference:
 Kerpedjiev *et al.*, 2018

Tool:
 JBrowse2

Version:
 {{ assembly.JBrowse2_version }}

Source:
 [https://jbrowse.org/jb2/](https://jbrowse.org/jb2/)

Reference:
 Diesh *et al.*, 2023

Tool:
 Juicer

Version:
 {{ assembly.Juicer_version }}

Source:
 [https://github.com/aidenlab/juicer](https://github.com/aidenlab/juicer)

Reference:
 Durand *et al.*, 2016

Tool:
 JuicerTools

Version:
 {{ assembly.JuicerTools_version }}

Source:
 [https://github.com/aidenlab/JuicerTools](https://github.com/aidenlab/JuicerTools)

Reference:
 Durand *et al.*, 2016 

Tool:
 Merqury.FK

Version:
 {{ assembly.MercuryFK_version }}

Source:
 [https://github.com/thegenemyers/MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

Reference:
 Rhie *et al.*, 2020

Tool:
 Minimap2

Version:
 {{ assembly.Minimap2_version }}

Source:
 [https://github.com/RemiAllio/MitoFinder](https://github.com/RemiAllio/MitoFinder)

Reference:
 Li, 2018

Tool:
 MitoFinder

Version:
 {{ assembly.MitoFinder_version }}

Source:
 [https://github.com/RemiAllio/MitoFinder](https://github.com/RemiAllio/MitoFinder)

Reference:
 Allio *et al.*, 2020

Tool:
 MitoHiFi

Version:
 {{ assembly.MitoHiFi_version }}

Source:
 [https://github.com/marcelauliano/MitoHiFi](https://github.com/marcelauliano/MitoHiFi)

Reference:
 Uliano-Silva *et al.*, 2023

Tool:
 MultiQC

Version:
 {{ assembly.MultiQC_version }}

Source:
 [https://github.com/MultiQC/MultiQC](https://github.com/MultiQC/MultiQC)

Reference:
 Ewels *et al.* 2016

Tool:
 NCBI Datasets

Version:
 {{ assembly.NCBI_Datasets_version }}

Source:
 [https://github.com/ncbi/datasets](https://github.com/ncbi/datasets)

Reference:
 Sayers *et al.*, 2024

Tool:
 Nextflow

Version:
 {{ assembly.Nextflow_version }}

Source:
 [https://github.com/nextflow-io/nextflow](https://github.com/nextflow-io/nextflow)

Reference:
 Di Tommaso *et al.*, 2017

Tool:
 Nf-core

Version:
 NA

Source:
 [https://nf-co.re/](https://nf-co.re/)

Reference:
 Ewels *et al.*, 2020

Tool:
 Oatk

Version:
 {{ assembly.oatk_version }}

Source:
 [https://github.com/c-zhou/oatk](https://github.com/c-zhou/oatk)

Reference:
 Zhou *et al.*, 2024

Tool:
 pigz

Version:
 {{ assembly.pigz_version }}

Source:
 [https://github.com/madler/pigz](https://github.com/madler/pigz)

Reference:
 -

Tool:
 PretextMap

Version:
 {{ assembly.PretextMap_version }}

Source:
 [https://github.com/sanger-tol/PretextMap](https://github.com/sanger-tol/PretextMap)

Reference:
 -

Tool:
 PretextSnapshot

Version:
 {{ assembly.PretextSnapshot_version }}

Source:
 [https://github.com/sanger-tol/PretextSnapshot](https://github.com/sanger-tol/PretextSnapshot)

Reference:
 -

Tool:
 purge_dups

Version:
 {{ assembly.purge_dups_version }}

Source:
 [https://github.com/dfguan/purge_dups](https://github.com/dfguan/purge_dups)

Reference:
 Guan *et al.*, 2020

Tool:
 samtools

Version:
 {{ assembly.samtools_version }}

Source:
 [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

Reference:
 Danecek *et al.*, 2021

Tool:
Seqtk

Version:
 {{ assembly.Seqtk_version }}

Source:
 [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

Reference:
 NA

Tool:
 Singularity

Version:
 {{ assembly.Singluarity_version }}

Source:
 [https://github.com/sylabs/singularity](https://github.com/sylabs/singularity)

Reference:
 Kurtzer *et al.*, 2017

Tool:
 Uniprot

Version:
 {{ assembly.Uniprot_version }}

Source:
 [https://www.uniprot.org/](https://www.uniprot.org/)

Reference:
 Bateman *et al.*, 2023

Tool:
 YaHS

Version:
 {{ assembly.YaHS_version }}

Source:
 [https://github.com/c-zhou/yahs](https://github.com/c-zhou/yahs)

Reference:
 Zhou *et al.*, 2023

## **Data availability**

Raw sequencing data, sample metadata, genome assembly sequences and
annotation data are available from the European Nucleotide Archive under
the BioProject accession number {{ bioproject_accession }};
[https://identifiers.org/ena.embl/](https://identifiers.org/ena.embl/){{
bioproject_accession }}. Assembly and raw data accession identifiers are
reported in Tables 1 and 4. Raw sequence data and sample metadata were
originally submitted to the Bioplatforms Australia Data Portal
([https://data.bioplatforms.com/](https://data.bioplatforms.com/)),
and are available under the following data package identifiers: {{
bpa_package_id }}.

The genome sequence is released openly for reuse.

## **Grant information**

Samples were collected and sequence data generated as part of the
Bioplatforms Australia-sponsored sequencing project {{ sample.bpa_initiative
}}. The genome has been assembled, annotated, and published as part of
the Australian Tree of Life Informatics Capacity, provided by the
Australian BioCommons.

## **Acknowledgements**

The pipelines used to assemble and publish this genome sequence have been adapted from the original digital infrastructure developed as part of the Darwin Tree of Life project (Darwin Tree of Life Project Consortium, 2022). This generation of this sequence report has leveraged assets from the Tree of Life Genome Note pipeline (Babirye *et al.*, 2025).

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
