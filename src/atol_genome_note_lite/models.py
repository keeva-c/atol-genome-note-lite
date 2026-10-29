from __future__ import annotations

import re
import sys
from datetime import date, datetime, time
from decimal import Decimal
from enum import Enum
from typing import Any, ClassVar, Literal, Optional, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    field_validator,
    model_serializer,
)

metamodel_version = "1.12.0"
version = "None"


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        serialize_by_alias=True,
        validate_by_name=True,
        validate_assignment=True,
        validate_default=True,
        extra="ignore",
        arbitrary_types_allowed=True,
        use_enum_values=True,
        strict=False,
    )


class LinkMLMeta(RootModel):
    root: dict[str, Any] = {}
    model_config = ConfigDict(frozen=True)

    def __getattr__(self, key: str):
        return getattr(self.root, key)

    def __getitem__(self, key: str):
        return self.root[key]

    def __setitem__(self, key: str, value):
        self.root[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self.root


linkml_meta = LinkMLMeta(
    {
        "default_prefix": "https://example.org/atol-report-schema/",
        "default_range": "string",
        "description": "this is a test of modelling the atol genome report schema "
        "with linkml",
        "id": "https://example.org/atol-report-schema",
        "imports": ["linkml:types", "busco_json_schema"],
        "name": "atol-genome-report-schema",
        "prefixes": {
            "linkml": {
                "prefix_prefix": "linkml",
                "prefix_reference": "https://w3id.org/linkml/",
            }
        },
        "source_file": "src/atol_genome_note_lite/schema/atol_genome_report_schema.yaml",
    }
)


class BuscoStats(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {
            "from_schema": "https://example.org/atol-report-schema-busco",
            "tree_root": True,
        }
    )

    parameters: Optional[BuscoParameter] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoStats"]}}
    )
    lineage_dataset: Optional[BuscoLineageDataset] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {"domain_of": ["BuscoStats", "BuscoParameter"]}
        },
    )
    versions: Optional[BuscoVersion] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoStats"]}}
    )
    results: Optional[BuscoResult] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoStats"]}}
    )
    metrics: Optional[BuscoMetric] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoStats"]}}
    )


class BuscoParameter(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema-busco"}
    )

    out_path: Optional[str] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    cpu: Optional[int] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    force: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    restart: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    quiet: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    download_path: Optional[str] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    datasets_version: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "odb12"}],
            }
        },
    )
    offline: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    download_base_url: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "https://busco-data.ezlab.org/v5/data/"}],
            }
        },
    )
    auto_lineage: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    auto_lineage_prok: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    auto_lineage_euk: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    use_augustus: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    use_miniprot: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    use_metaeuk: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    batch_mode: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    tar: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    opt_out_run_stats: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    plot_percentages: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    out: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "asm-busco"}],
            }
        },
    )
    mode: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "euk_genome_min"}],
            }
        },
    )
    lineage_dataset: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoStats", "BuscoParameter"],
                "examples": [{"value": "/path/to/lineages/lepidoptera_odb12"}],
            }
        },
    )
    ncbi_taxid: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "7088"}],
            }
        },
    )
    orthodb_version: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "12.1"}],
            }
        },
    )
    dataset_version: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "02"}],
            }
        },
    )
    contig_break: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "10"}],
            }
        },
    )
    scaffold_composition: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    in_: Optional[str] = Field(
        default=None,
        alias="in",
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "/path/to/asm_hap1_scaffolds_final.fa"}],
            }
        },
    )
    main_out: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "/path/to/asm_hap1_scaffolds_final.fa"}],
            }
        },
    )
    lineage_results_dir: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "run_lepidoptera_odb12"}],
            }
        },
    )
    domain: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter", "BuscoResult"],
                "examples": [{"value": "eukaryota"}],
            }
        },
    )
    skip_bbtools: Optional[bool] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["BuscoParameter"]}},
    )
    gene_predictor: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter"],
                "examples": [{"value": "miniprot"}],
            }
        },
    )


class BuscoLineageDataset(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema-busco"}
    )

    name: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoLineageDataset"],
                "examples": [{"value": "lepidoptera_odb12"}],
            }
        },
    )
    creation_date: Optional[date] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoLineageDataset"],
                "examples": [{"value": "2026-05-22"}],
            }
        },
    )
    number_of_buscos: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoLineageDataset"],
                "examples": [{"value": "5760"}],
            }
        },
    )
    number_of_species: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoLineageDataset"],
                "examples": [{"value": "79"}],
            }
        },
    )


class BuscoVersion(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema-busco"}
    )

    python: Optional[list[Union[int, str]]] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "any_of": [{"range": "integer"}, {"range": "string"}],
                "domain_of": ["BuscoVersion"],
                "examples": [{"value": "[3, 12, 11, 'final', 0]"}],
            }
        },
    )
    busco: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoVersion"],
                "examples": [{"value": "6.0.0"}],
            }
        },
    )
    hmmsearch: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoVersion"],
                "examples": [{"value": "3.4"}],
            }
        },
    )
    bbtools: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoVersion"]}}
    )
    miniprot_index: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoVersion"],
                "examples": [{"value": "0.17-r279"}],
            }
        },
    )
    miniprot_align: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoVersion"],
                "examples": [{"value": "0.17-r279"}],
            }
        },
    )


class BuscoResult(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema-busco"}
    )

    one_line_summary: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [
                    {"value": "C:99.0%[S:96.6%,D:2.4%],F:0.1%,M:0.9%,n:5760,E:7.3%"}
                ],
            }
        },
    )
    Complete_percentage: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "99.0"}],
            }
        },
    )
    Complete_BUSCOs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "5703"}],
            }
        },
    )
    Single_copy_percentage: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "96.6"}],
            }
        },
    )
    Single_copy_BUSCOs: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["BuscoResult"]}}
    )
    Multi_copy_percentage: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "2.4"}],
            }
        },
    )
    Multi_copy_BUSCOs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "139"}],
            }
        },
    )
    Fragmented_percentage: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "0.1"}],
            }
        },
    )
    Fragmented_BUSCOs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {"domain_of": ["BuscoResult"], "examples": [{"value": "5"}]}
        },
    )
    Missing_percentage: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "0.9"}],
            }
        },
    )
    Missing_BUSCOs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {"domain_of": ["BuscoResult"], "examples": [{"value": "52"}]}
        },
    )
    n_markers: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "5760"}],
            }
        },
    )
    domain: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoParameter", "BuscoResult"],
                "examples": [{"value": "eukaryota"}],
            }
        },
    )
    avg_identity: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "0.81"}],
            }
        },
    )
    internal_stop_codon_count: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "419"}],
            }
        },
    )
    internal_stop_codon_percent: Optional[float] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult"],
                "examples": [{"value": "7.3"}],
            }
        },
    )
    Number_of_scaffolds: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "447"}],
            }
        },
    )
    Number_of_contigs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "465"}],
            }
        },
    )
    Total_length: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "531458053"}],
            }
        },
    )
    Percent_gaps: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "0.000%"}],
            }
        },
    )
    Scaffold_N50: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "17662296"}],
            }
        },
    )
    Contigs_N50: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "14578451"}],
            }
        },
    )


class BuscoMetric(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema-busco"}
    )

    Number_of_scaffolds: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "447"}],
            }
        },
    )
    Number_of_contigs: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "465"}],
            }
        },
    )
    Total_length: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "531458053"}],
            }
        },
    )
    Percent_gaps: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "0.000%"}],
            }
        },
    )
    Scaffold_N50: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "17662296"}],
            }
        },
    )
    Contigs_N50: Optional[int] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["BuscoResult", "BuscoMetric"],
                "examples": [{"value": "14578451"}],
            }
        },
    )


class TaxonomyInfo(ConfiguredBaseModel):
    """
    Containing information retrieved from Canopy's Taxonomy table
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    busco_odb12_dataset_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_authority: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_common_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_family: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_full_lineage: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_order: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_scientific_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )
    ncbi_taxon_id: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["TaxonomyInfo"]}}
    )


class SeqStats(ConfiguredBaseModel):
    """
    A line of results from seqkit stats
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    avg_len: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    AvgQual: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    format: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    GC___: Optional[float] = Field(
        default=None,
        alias="GC(%)",
        json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}},
    )
    max_len: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    min_len: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    N50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    N50_num: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    num_seqs: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    Q1: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    Q2: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    Q20___: Optional[int] = Field(
        default=None,
        alias="Q20(%)",
        json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}},
    )
    Q3: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    Q30___: Optional[int] = Field(
        default=None,
        alias="Q30(%)",
        json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}},
    )
    sequencing_depth: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    sum_gap: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    sum_len: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    sum_n: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    type: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["SeqStats"]}}
    )
    sum_len_print: Optional[str] = Field(
        default=None,
        json_schema_extra={
            "linkml_meta": {
                "domain_of": ["SeqStats"],
                "string_serialization": "{sum_len} // 1000",
            }
        },
    )


class MerquryStats(ConfiguredBaseModel):
    """
    A line of QV and Completeness stats from Merqury
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    Region: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    Found: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    Total: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    Covered: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    No_Support: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    Error: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )
    QV: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["MerquryStats"]}}
    )


class SoftwareVersions(ConfiguredBaseModel):
    """
    The versions of pipelines used in assembly generation
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    sanger_tol_genomeassembly: Optional[str] = Field(
        default=None,
        alias="sanger_tol/genomeassembly",
        json_schema_extra={"linkml_meta": {"domain_of": ["SoftwareVersions"]}},
    )
    sanger_tol_ascc: Optional[str] = Field(
        default=None,
        alias="sanger_tol/ascc",
        json_schema_extra={"linkml_meta": {"domain_of": ["SoftwareVersions"]}},
    )
    sanger_tol_treeval: Optional[str] = Field(
        default=None,
        alias="sanger_tol/treeval",
        json_schema_extra={"linkml_meta": {"domain_of": ["SoftwareVersions"]}},
    )


class PackageDataContainer(ConfiguredBaseModel):
    """
    Container class containing all metadata and metrics associated with one data package on the Bioplatforms data portal
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema", "tree_root": True}
    )

    project: Optional[Project] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    organism: Optional[Organism] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    specimen: Optional[Specimen] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    sample: Optional[Sample] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    experiment: Optional[Experiment] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    runs: Optional[Runs] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    assembly: Optional[Assembly] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )
    annotation: Optional[Annotation] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["PackageDataContainer"]}},
    )


class Project(ConfiguredBaseModel):
    """
    Containing information related to the 'root' level bioproject
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    project_accession: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Project"]}}
    )


class Organism(ConfiguredBaseModel):
    """
    Containing information related to the taxon
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    taxonomy_info: Optional[TaxonomyInfo] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Organism"]}}
    )


class Specimen(ConfiguredBaseModel):
    """
    Containing information related to the specimen-level sample
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    biosample_accession: Optional[str] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["Specimen", "Sample"]}},
    )


class Sample(ConfiguredBaseModel):
    """
    Containing information related to the sample-level sample
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    project_collaborators: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    bpa_sample_id: str = Field(
        default=..., json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    bioplatforms_project_id: str = Field(
        default=..., json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    bpa_initiative: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    biosample_accession: Optional[str] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["Specimen", "Sample"]}},
    )
    specimen_id: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    specimen_id_description: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    sex: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    lifestage: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    collection_date: Optional[date] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    region_and_locality: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    state_or_region: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    country_or_sea: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    indigenous_location: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    latitude: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    longitude: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    collected_by: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    collection_method: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    collection_permit: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    identified_by: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    preservation_method: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    preservation_temperature: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )
    organism_part: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Sample"]}}
    )


class Experiment(ConfiguredBaseModel):
    """
    Containing information related to the experiment (ERX)
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    data_owner: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    bpa_library_id: str = Field(
        default=..., json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    library_strategy: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    extraction_method: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    nucleic_acid_treatment: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    nucleic_acid_conc: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    library_construction_protocol: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    library_selection: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    insert_size: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    platform: str = Field(
        default=..., json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    instrument_model: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    sequencing_kit: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    gal: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    library_layout: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    flowcell_type: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    base_caller_model: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )
    bpa_package_id: str = Field(
        default=..., json_schema_extra={"linkml_meta": {"domain_of": ["Experiment"]}}
    )


class Runs(ConfiguredBaseModel):
    """
    Containing information related to the qc reads (ERR)
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    submission_records: Optional[SubmissionRecords] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Runs"]}}
    )
    read_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Runs"]}}
    )
    base_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Runs"]}}
    )


class Assembly(ConfiguredBaseModel):
    """
    Containing information related to the genome assembly
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    assembly_level: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    haplotypes: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    rna_data_available: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    tol_id: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    version: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    assembly_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_assembly_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_assembly_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    assembly_accession: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_assembly_accession: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_assembly_accession: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    contaminants_removed: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_contaminants_removed: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_contaminants_removed: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    coverage: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    scaffold_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    contig_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    genome_length: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    scaffold_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    contig_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    gap_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    longest_contig: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    longest_scaffold: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    busco_string: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    busco_single: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    busco_duplicated: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    primary_qv: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    alt_qv: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    combined_qv: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    primary_kmer: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    alt_kmer: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    combined_kmer: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    mito_size: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    plastid_size: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_scaffold_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_genome_length: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_scaffold_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_contig_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_busco_c: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_gap_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_contig_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_longest_contig: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_longest_scaffold: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_busco_string: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_busco_single: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_busco_duplicated: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_scaffold_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_genome_length: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_scaffold_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_contig_n50: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_busco_c: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_gap_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_contig_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_longest_contig: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_longest_scaffold: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_busco_string: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_busco_single: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_busco_duplicated: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    genomescope_image_path: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    contact_map_image_path: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_1_contact_map_image_path: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    hap_2_contact_map_image_path: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    genomeassembly_pipeline_version: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    ascc_pipeline_version: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )
    treeval_pipeline_version: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Assembly"]}}
    )


class Annotation(ConfiguredBaseModel):
    """
    Containing information related to the genome Annotation
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    annotation_pipeline: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annotation_pipeline_ver: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annotation_pipeline_link: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annot_qc_pipeline: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annot_qc_pipeline_ver: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annot_qc_pipeline_link: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annotation_name: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    gene_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    cds_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    transcript_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    mean_transcript_length: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    mean_transcripts_per_gene: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    mean_exons_per_transcript: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annot_busco_summary: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    annot_busco_lineage: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_completeness_summary: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    conserved_hogs: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_lineage: Optional[str] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_protein_count: Optional[int] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_percent_consistent: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_percent_inconsistent: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_percent_contaminant: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )
    omark_percent_unknown: Optional[float] = Field(
        default=None, json_schema_extra={"linkml_meta": {"domain_of": ["Annotation"]}}
    )


class SubmissionRecords(ConfiguredBaseModel):
    """
    Containing information about the ERR submission record
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {"from_schema": "https://example.org/atol-report-schema"}
    )

    accession: Optional[str] = Field(
        default=None,
        json_schema_extra={"linkml_meta": {"domain_of": ["SubmissionRecords"]}},
    )


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
BuscoStats.model_rebuild()
BuscoParameter.model_rebuild()
BuscoLineageDataset.model_rebuild()
BuscoVersion.model_rebuild()
BuscoResult.model_rebuild()
BuscoMetric.model_rebuild()
TaxonomyInfo.model_rebuild()
SeqStats.model_rebuild()
MerquryStats.model_rebuild()
SoftwareVersions.model_rebuild()
PackageDataContainer.model_rebuild()
Project.model_rebuild()
Organism.model_rebuild()
Specimen.model_rebuild()
Sample.model_rebuild()
Experiment.model_rebuild()
Runs.model_rebuild()
Assembly.model_rebuild()
Annotation.model_rebuild()
SubmissionRecords.model_rebuild()
