#!/usr/bin/env python3

import argparse
from importlib.resources import files
from pathlib import Path

from black import FileMode, format_str
from linkml.generators.pydanticgen import PydanticGenerator


def _existing_file(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return p


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

    models = PydanticGenerator(
        schema=args.schema,
        black=True,
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
