#!/usr/bin/env python3

import logging
from pathlib import Path


def _existing_file(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return p


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

logger = logging.getLogger()
