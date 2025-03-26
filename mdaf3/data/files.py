"""
Location of data files
======================

Use as ::

    from mdaf3.data.files import *

"""

__all__ = [
    "UNCOMPRESSED_AF3_OUTPUT_PATH",
]

import importlib.resources

data_directory = importlib.resources.files("mdaf3") / "data"

UNCOMPRESSED_AF3_OUTPUT_PATH = (
    data_directory / "93f0240a1d2c15da9551841d22239d41"
)
