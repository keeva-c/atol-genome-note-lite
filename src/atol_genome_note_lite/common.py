#!/usr/bin/env python3


def _existing_file(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return p
