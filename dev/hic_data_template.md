# **extra sample and experiment info to add to the genome note (Hi-C)**

**Sample information**

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


**Collection information**  

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


**Identification information**

Data type generated\
 Hi-C

Identified by\
 {{ sample.identified_by }}

**Preservation information**

Data type generated\
 Hi-C

Preservation method\
 {{ sample.preservation_method }}

Preservation temperature\
 {{ sample.preservation_temperature }}


## **Nucleic acid extraction**

Table 3: Methodological information about nucleic acid material
extracted for sequencing.

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

## **Sequencing**

Table 4: Methodological information about sequencing runs.

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
