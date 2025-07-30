"""
Location of data files
======================

Use as ::

    from mdaf3.data.files import *

"""

__all__ = [
    "UNCOMPRESSED_AF3_OUTPUT_PATH",
    "UNCOMPRESSED_BOLTZ_OUTPUT_PATH",
]

import importlib.resources

data_directory = importlib.resources.files("mdaf3.data")

UNCOMPRESSED_AF3_OUTPUT_PATH = (
    data_directory / "4af4e6b7ce5cfda8793e03ba174b6be4"
)

UNCOMPRESSED_BOLTZ_OUTPUT_PATH = (
    data_directory / "1d3df2ab5ff6b4050b4a4abae1553a6d"
)
