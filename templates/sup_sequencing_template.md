| **Run: {{ runs.submission_records.accession }}** | |
| Data type generated | {{ experiment.library_strategy }} |
| Sample-level BioSample accession | {{ sample.biosample_accession }} |
| Library identifier | {{ experiment.bpa_library_id }} |
| Run accession | {{ runs.submission_records.accession }} |
| Read count | {% if runs.read_count %}{{ make_pretty_number(runs.read_count) }}{% else %}*not provided*{% endif %} |
| Base count | {% if runs.base_count %}{{ round_bases_up(runs.base_count) }}{% else %}*not provided*{% endif %} |
| Sequencing platform | {{ experiment.platform }} |
| Sequencing instrument | {{ experiment.instrument_model }} |
| Sequencing chemistry | {{ experiment.sequencing_kit }} |
| Sequencing facility | {{ experiment.gal }} |
| Library layout | {{ experiment.library_layout }} |
| Flowcell type | {{ experiment.flowcell_type }} | {% if experiment.platform=='Oxford Nanopore' %}
| Base caller model | {{ experiment.base_caller_model }} | {% endif %}