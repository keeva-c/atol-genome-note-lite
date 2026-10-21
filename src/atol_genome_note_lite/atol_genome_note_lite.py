#!/usr/bin/env python3

import argparse
from importlib.resources import files
import json
from pathlib import Path

from jinja2 import Environment

from atol_genome_note_lite.common import _existing_file
from atol_genome_note_lite.models import TaxonomyInfo
from canopy_client import CanopySession


def parse_arguments() -> argparse.Namespace:
    my_files = files()

    parser = argparse.ArgumentParser(
        description="Render the templates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    _ = parser.add_argument(
        "--title_template",
        type=_existing_file,
        default=Path(my_files, "templates", "title.md.j2"),
        help="Title template.",
    )

    return parser.parse_args()


def render_template_string(template_string: str, **kwargs) -> str:
    env = Environment()
    template = env.from_string(template_string)
    context = {**kwargs}
    return template.render(context)


def render_template(template_path: Path, **kwargs) -> str:
    return render_template_string(template_path.read_text(), **kwargs)


def main():

    args = parse_arguments()

    canopy_session = CanopySession()
    taxonomy_info = canopy_session.get_taxonomy_info(taxon_id=319761).json()

    # raise ValueError(json.dumps(taxonomy_info, indent=4))

    x = TaxonomyInfo(**taxonomy_info)

    raise ValueError(
        render_template(args.title_template, **{"taxonomy_info": x.model_dump()})
    )


if __name__ == "__main__":
    main()
