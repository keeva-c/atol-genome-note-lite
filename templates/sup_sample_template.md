| **Sample: {{ sample.biosample_accession }}** | |
| BioSample accession | {{ sample.biosample_accession }} |
| Scientific name | *{{ taxonomy_info.ncbi_scientific_name }}* |
| ToLID | {{ sample.tolid }} |
| Specimen identifier | {{ sample.specimen_id }} |
| Specimen identifier defined by | {{ sample.specimen_id_description }} |
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