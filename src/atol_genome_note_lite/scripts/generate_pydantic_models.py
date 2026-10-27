#!/usr/bin/env python3

import argparse
from importlib.resources import files
from pathlib import Path

from atol_genome_note_lite.common import _existing_file
from black import FileMode, format_str
from linkml.generators.pydanticgen import PydanticGenerator


def parse_arguments() -> argparse.Namespace:
    my_files = files()

    parser = argparse.ArgumentParser(
        description="Generate pydantic models",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    _ = parser.add_argument(
        "--schema",
        type=_existing_file,
        default=Path(my_files.parent, "schema", "atol_genome_report_schema.yaml"),
        help="Schema input.",
    )
    _ = parser.add_argument(
        "--models",
        type=Path,
        help="Pydantic models output.",
        default=Path(my_files.parent, "models.py"),
    )

    _ = parser.add_argument(
        "--force",
        help="Overwrite the Pydantic models output if it exists.",
        action="store_true",
    )

    _ = parser.add_argument("-n", "--dry_run", action="store_true", default=False)

    return parser.parse_args()


def main():

    args = parse_arguments()

    # NOTE. The PydanticGenerator generates aliases for slots in the schema
    # that contain special characters. e.g. the "GC(%)" slot from seqkit will
    # be called `GC___` with `alias="GC(%)"`. However, this doesn't seem to
    # handle whitespace, so "Complete BUSCOs" appears in the models as
    # `Complete_BUSCOs` with no alias. Because the BUSCO JSON uses `Complete
    # BUSCOs` as a key, we have to handle the whitespace when we load the data.
    # See the read_busco_json function in
    # src/atol_genome_note_lite/atol_genome_note_lite.py for an example.
    models = PydanticGenerator(
        schema=args.schema, black=True, use_aliases=True, extra_fields="ignore", 
    )

    formatted = format_str(models.serialize(), mode=FileMode())

    if args.dry_run is True:
        print(formatted)

    if args.models.exists() and args.force is False and args.dry_run is False:
        raise FileExistsError(args.models)

    if args.dry_run is False:
        args.models.parent.mkdir(parents=True, exist_ok=True)
        with open(args.models, "wt") as f:
            _ = f.write(formatted)


if __name__ == "__main__":
    main()
