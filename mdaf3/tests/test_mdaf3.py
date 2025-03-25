"""
Unit and regression test for the mdaf3 package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import mdaf3


def test_mdaf3_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mdaf3" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"
