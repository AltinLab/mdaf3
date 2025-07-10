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


class TestAF3CompressedUnchanged:

    @pytest.fixture
    def uncompressed_output(self):
        yield AF3Output(UNCOMPRESSED_AF3_OUTPUT_PATH)

    @pytest.fixture
    def compressed_output(self, tmp_path):
        compressed_af3_output_path = (
            tmp_path / UNCOMPRESSED_AF3_OUTPUT_PATH.name
        )
        copyanything(UNCOMPRESSED_AF3_OUTPUT_PATH, compressed_af3_output_path)

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

        # AF3 only stores 2 decimal points
        assert_allclose(comp_pae_ndarr, uncomp_pae_ndarr, atol=0.01)

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

        assert_allclose(
            comp_contact_prob_ndarr, uncomp_contact_prob_ndarr, atol=0.01
        )

    @pytest.mark.parametrize("seed", [None, 1, 2])
    def test_summary_metrics(
        self, uncompressed_output, compressed_output, seed
    ):
        comp_summary_metrics = compressed_output.get_summary_metrics(seed=seed)
        uncomp_summary_metrics = uncompressed_output.get_summary_metrics(
            seed=seed
        )

        for k in comp_summary_metrics.keys():
            assert np.all(
                np.allclose(
                    comp_summary_metrics[k],
                    uncomp_summary_metrics[k],
                    atol=0.01,
                )
            )

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
        comp_res_ids = compressed_output.get_token_res_ids(seed=seed)
        uncomp_res_ids = uncompressed_output.get_token_res_ids(seed=seed)

        assert np.all(np.equal(comp_res_ids, uncomp_res_ids))


@pytest.mark.parametrize("seed", [None, 1, 2])
@pytest.mark.parametrize("compressed", [False, True])
class TestFeatureSelection:

    @pytest.fixture
    def uncompressed_output(self):
        yield AF3Output(UNCOMPRESSED_AF3_OUTPUT_PATH)

    @pytest.fixture
    def compressed_output(self, tmp_path):
        compressed_af3_output_path = (
            tmp_path / UNCOMPRESSED_AF3_OUTPUT_PATH.name
        )
        copyanything(UNCOMPRESSED_AF3_OUTPUT_PATH, compressed_af3_output_path)

        compressed_af3_output = AF3Output(compressed_af3_output_path)
        compressed_af3_output.compress()

        if (
            compressed_af3_output.compressed == False
            or (compressed_af3_output.dir_path / "TERMS_OF_USE.md").is_file()
        ):
            raise ValueError("Compression failed")

        yield compressed_af3_output

    @pytest.fixture
    def af3_output(self, compressed, compressed_output, uncompressed_output):
        if compressed:
            yield compressed_output
        else:
            yield uncompressed_output

    @pytest.fixture
    def universe(self, af3_output, seed):
        yield af3_output.get_mda_universe(seed=seed)

    def test_chain_pair_pae_min(self, universe, af3_output, seed):
        chain_pair_pae_min = af3_output.get_summary_metrics(seed=seed)[
            "chain_pair_pae_min"
        ]

        pae = af3_output.get_pae_ndarr(seed=seed)
        chain_i_resindices = []
        for segid in ["A", "B", "C", "D", "E"]:
            chain_i_resindices.append(
                universe.select_atoms(f"segid {segid}").residues.resindices
            )

        manual_feat = np.zeros((5, 5), dtype=np.float32)
        for i in range(5):
            for j in range(5):
                manual_feat[i][j] = pae[chain_i_resindices[i]][
                    :, chain_i_resindices[j]
                ].min()
                manual_feat[j][i] = pae[chain_i_resindices[j]][
                    :, chain_i_resindices[i]
                ].min()
        # since PAE is only stored to 2 decimal points,
        # but chain_pair_pae_min is calculated pre-truncation
        # we can't assert to a high tolerance here
        assert_allclose(chain_pair_pae_min, manual_feat, atol=0.1)
