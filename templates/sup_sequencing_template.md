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