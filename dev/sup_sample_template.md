**Sample information**

BioSample accession\
 {{ sample.biosample_accession }}

Scientific name\
 *{{ organism.scientific_name }}* 

TOLID\
 {{ organism.tolid }}

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