"""
Unit and regression test for the mdaf3 package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import mdaf3
import shutil, errno
from mdaf3.data.files import *
from mdaf3.AF3OutputParser import *
from numpy.testing import assert_allclose


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno in (errno.ENOTDIR, errno.EINVAL):
            shutil.copy(src, dst)
        else:
            raise


class TestCompressedUnchanged:

    @pytest.fixture(scope="class")
    def uncompressed_output(self):
        yield AF3Output(UNCOMPRESSEED_AF3_OUTPUT_PATH)

    @pytest.fixture(scope="class")
    def compressed_output(self, tmp_path):
        compressed_af3_output_path = (
            tmp_path / UNCOMPRESSEED_AF3_OUTPUT_PATH.name
        )
        copyanything(UNCOMPRESSEED_AF3_OUTPUT_PATH, compressed_af3_output_path)

        compressed_af3_output = AF3Output(compressed_af3_output_path)
        compressed_af3_output.compress()

        if (
            compressed_af3_output.compressed == False
            or (compressed_af3_output.dir_path / "TERMS_OF_USE.md").is_file()
        ):
            raise ValueError("Compression failed")

        yield compressed_af3_output

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_pae_ndarr(self, uncompressed_output, compressed_output, seed):
        comp_pae_ndarr = compressed_output.get_pae_ndarr(seed=seed)
        uncomp_pae_ndarr = uncompressed_output.get_pae_ndarr(seed=seed)

        assert_allclose(comp_pae_ndarr, uncomp_pae_ndarr)

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_contact_prob_ndarr(
        self, uncompressed_output, compressed_output, seed
    ):
        comp_contact_prob_ndarr = compressed_output.get_contact_prob_ndarr(
            seed=seed
        )
        uncomp_contact_prob_ndarr = uncompressed_output.get_contact_prob_ndarr(
            seed=seed
        )

        assert_allclose(comp_contact_prob_ndarr, uncomp_contact_prob_ndarr)

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_summary_metrics(
        self, uncompressed_output, compressed_output, seed
    ):
        comp_summary_metrics = compressed_output.get_summary_metrics(seed=seed)
        uncomp_summary_metrics = uncompressed_output.get_summary_metrics(
            seed=seed
        )

        assert comp_summary_metrics == uncomp_summary_metrics

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_token_chain_ids(
        self, uncompressed_output, compressed_output, seed
    ):
        comp_chain_ids = compressed_output.get_token_chain_ids(seed=seed)
        uncomp_chain_ids = uncompressed_output.get_token_chain_ids(seed=seed)

        assert np.all(np.equal(comp_chain_ids, uncomp_chain_ids))

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_atom_chain_ids(
        self, uncompressed_output, compressed_output, seed
    ):
        comp_chain_ids = compressed_output.get_atom_chain_ids(seed=seed)
        uncomp_chain_ids = uncompressed_output.get_atom_chain_ids(seed=seed)

        assert np.all(np.equal(comp_chain_ids, uncomp_chain_ids))

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_token_res_ids(self, uncompressed_output, compressed_output, seed):
        comp_res_ids = compressed_output.get_res_ids(seed=seed)
        uncomp_res_ids = uncompressed_output.get_res_ids(seed=seed)

        assert np.all(np.equal(comp_res_ids, uncomp_res_ids))


@pytest.mark.parametrize("seed", [None, 1, 2])
@pytest.mark.parametrize("compressed", [False, True])
class TestFeatureSelection:

    @pytest.fixture(scope="class")
    def uncompressed_output(self):
        yield AF3Output(UNCOMPRESSEED_AF3_OUTPUT_PATH)

    @pytest.fixture(scope="class")
    def compressed_output(self, tmp_path):
        compressed_af3_output_path = (
            tmp_path / UNCOMPRESSEED_AF3_OUTPUT_PATH.name
        )
        copyanything(UNCOMPRESSEED_AF3_OUTPUT_PATH, compressed_af3_output_path)

        compressed_af3_output = AF3Output(compressed_af3_output_path)
        compressed_af3_output.compress()

        if (
            compressed_af3_output.compressed == False
            or (compressed_af3_output.dir_path / "TERMS_OF_USE.md").is_file()
        ):
            raise ValueError("Compression failed")

        yield compressed_af3_output

    @pytest.fixture
    def af3_output(self, compressed):
        if compressed:
            yield compressed_output
        else:
            yield uncompressed_output
