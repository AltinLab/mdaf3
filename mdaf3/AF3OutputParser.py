import numpy as np
import orjson
from pathlib import Path
import polars as pl
import MDAnalysis as mda
import polars as pl
import h5py
import gzip

SUMMARY_ARRAY = [
    "chain_iptm",
    "chain_pair_iptm",
    "chain_pair_pae_min",
    "chain_ptm",
]

SUMMARY_SCALAR = [
    "fraction_disordered",
    "has_clash",
    "iptm",
    "ptm",
    "ranking_score",
]


class AF3Output:
    """
    Handle for the output directory of an AlphaFold3 job.
    """

    def __init__(self, directory_path, server=False, **kwargs):

        dir_path = None
        if isinstance(directory_path, str):
            dir_path = Path(directory_path)
        elif isinstance(directory_path, Path):
            dir_path = directory_path
        else:
            raise ValueError

        self.job_name = dir_path.name
        self.dir_path = dir_path
        self.server = server

        self.compressed_file = Path(self.dir_path / f"{self.job_name}.h5")

        if self.compressed_file.is_file():
            self.compressed = True
        else:
            self.compressed = False

        if server:
            self._best_model_path = dir_path / f"{self.job_name}_model_0.cif"
            self._best_full_data_path = (
                dir_path / f"{self.job_name}_full_data_0.json"
            )
            self._best_summary_path = (
                dir_path / f"{self.job_name}_summary_confidences_0.json"
            )
            raise NotImplementedError

        else:
            self._best_model_path = dir_path / f"{self.job_name}_model.cif"
            self._best_full_data_path = (
                dir_path / f"{self.job_name}_confidences.json"
            )
            self._best_summary_path = (
                dir_path / f"{self.job_name}_summary_confidences.json"
            )

    def get_contact_prob_ndarr(
        self, seed=None, sample_num=None, indices=(slice(None), slice(None))
    ):
        """
        (n_tokens, n_tokens) ndarray of contact probabilities
        """

        if self.compressed:
            with self.get_h5_handle() as hf:
                contact_probs = hf[
                    self._seed_samplenum_str(seed, sample_num)
                    + "/"
                    + "confidences"
                ]["contact_probs"]

                return (contact_probs[indices[0]][:, indices[1]] / 100).astype(
                    np.float32
                )

        else:

            full_data_dict = self._get_full_data_dict(
                seed=seed, sample_num=sample_num
            )
            contact_probs = np.array(
                full_data_dict["contact_probs"], dtype=np.float32
            )

            return contact_probs[indices[0]][:, indices[1]]

    def get_pae_ndarr(
        self, seed=None, sample_num=None, indices=(slice(None), slice(None))
    ):
        """
        (n_tokens, n_tokens) ndarray of predicted aligment errors
        """

        if self.compressed:
            with self._get_h5_handle() as hf:
                pae = hf[
                    self._seed_samplenum_str(seed, sample_num)
                    + "/"
                    + "confidences"
                ]["pae"]

                return (pae[indices[0]][:, indices[1]] / 100).astype(
                    np.float32
                )

        else:
            full_data_dict = self._get_full_data_dict(
                seed=seed, sample_num=sample_num
            )
            pae = np.array(full_data_dict["pae"], dtype=np.float32)

            return pae[indices[0]][:, indices[1]]

    def get_summary_metrics(self, seed=None, sample_num=None):
        if self.compressed:
            with self.get_h5_handle() as hf:
                summary_metrics = hf[
                    self._seed_samplenum_str(seed, sample_num)
                    + "/"
                    + "summary_confidences"
                ]

                return {
                    k: summary_metrics[k][()]
                    for k in SUMMARY_ARRAY + SUMMARY_SCALAR
                }
        else:
            return {
                k: np.array(v)[()]
                for k, v in self._get_summary_data_dict(
                    seed=seed, sample_num=sample_num
                ).items()
            }

    def get_mda_universe(self, seed=None, sample_num=None, **kwargs):
        """
        Get MDAnalysis universe for a sample

        If seed and sample_num are None, return the best model
        """
        subdir = self._seed_samplenum_str(seed, sample_num)

        if subdir == "":
            top_path = self._best_model_path
        else:
            top_path = self.dir_path / subdir / "model.cif"

        if self.compressed:
            top_path = self._append_gz(top_path)

        return mda.Universe(
            top_path.as_posix(), topology_format="MMCIF", **kwargs
        )

    def compress(self):
        if self.server:
            raise NotImplementedError
        if self.compressed:
            return

        ranking_scores = pl.read_csv(self.dir_path / "ranking_scores.csv")

        seed_np = ranking_scores.select("seed").to_series().to_numpy()
        sample_np = ranking_scores.select("sample").to_series().to_numpy()
        ranking_score_np = (
            ranking_scores.select("ranking_score").to_series().to_numpy()
        )

        with h5py.File((self.dir_path / f"{self.job_name}.h5"), "w") as hf:

            rank_grp = hf.create_group("ranking_scores")

            rank_grp.create_dataset("seed", data=seed_np)
            rank_grp.create_dataset("sample", data=sample_np)
            rank_grp.create_dataset("ranking_score", data=ranking_score_np)

            for seed, sample in zip(seed_np, sample_np):
                subdir = self._seed_samplenum_str(seed, sample)

                full_data_file = self.dir_path / subdir / "confidences.json"
                full_data_dict = orjson.loads(full_data_file.read_text())
                summary_file = (
                    self.dir_path / subdir / "summary_confidences.json"
                )
                summary_dict = orjson.loads(summary_file.read_text())

                grp = hf.create_group(subdir)

                conf_grp = grp.create_group("confidences")
                summ_grp = grp.create_group("summary_confidences")

                conf_grp.create_dataset(
                    "atom_chain_ids",
                    data=full_data_dict["atom_chain_ids"],
                    dtype="S1",
                    compression="gzip",
                    compression_opts=9,
                )

                conf_grp.create_dataset(
                    "token_chain_ids",
                    data=full_data_dict["token_chain_ids"],
                    dtype="S1",
                    compression="gzip",
                    compression_opts=9,
                )

                conf_grp.create_dataset(
                    "token_res_ids",
                    data=full_data_dict["token_res_ids"],
                    dtype=np.int16,
                    compression="gzip",
                    compression_opts=9,
                )

                conf_grp.create_dataset(
                    "contact_probs",
                    data=np.array(full_data_dict["contact_probs"]) * 100,
                    dtype=np.int16,
                    compression="gzip",
                    compression_opts=9,
                )
                conf_grp.create_dataset(
                    "pae",
                    data=np.array(full_data_dict["pae"]) * 100,
                    dtype=np.int16,
                    compression="gzip",
                    compression_opts=9,
                )
                conf_grp.create_dataset(
                    "atom_plddts",
                    data=np.array(full_data_dict["atom_plddts"]) * 100,
                    dtype=np.int16,
                    compression="gzip",
                    compression_opts=9,
                )

                for summary_stat in SUMMARY_ARRAY:
                    summ_grp.create_dataset(
                        summary_stat,
                        data=summary_dict[summary_stat],
                        dtype=np.float32,
                        compression="gzip",
                        compression_opts=9,
                    )

                for summary_stat in SUMMARY_SCALAR:
                    summ_grp.create_dataset(
                        summary_stat,
                        data=summary_dict[summary_stat],
                        dtype=np.float32,
                    )

                # compress all CIFS
                cif_path = self.dir_path / subdir

                with (
                    open(cif_path / "model.cif", "rb") as src,
                    gzip.open(cif_path / "model.cif.gz", "wb") as dst,
                ):
                    dst.writelines(src)

                # remove uncompressed CIFS
                (cif_path / "model.cif").unlink(missing_ok=True)

                # remove json files
                full_data_file.unlink(missing_ok=True)
                summary_file.unlink(missing_ok=True)

            # best dataset
            max_idx = np.argmax(ranking_score_np)

            # hard link to the best dataset
            best_dataset_dir = hf[
                f"seed-{seed_np[max_idx]}_sample-{sample_np[max_idx]}"
            ]
            hf["confidences"] = best_dataset_dir["confidences"]
            hf["summary_confidences"] = best_dataset_dir["summary_confidences"]

        # compress best cif
        with (
            open(self._best_model_path, "rb") as src,
            gzip.open(
                self._append_gz(self._best_model_path),
                "wb",
            ) as dst,
        ):
            dst.writelines(src)

        # remove uncompressed best cif
        self._best_model_path.unlink(missing_ok=True)
        # remove best json
        self._best_full_data_path.unlink(missing_ok=True)
        self._best_summary_path.unlink(missing_ok=True)

        # remove ranking score
        (self.dir_path / "ranking_scores.csv").unlink(missing_ok=True)

        (self.dir_path / f"{self.job_name}_data.json").unlink(missing_ok=True)
        (self.dir_path / "TERMS_OF_USE.md").unlink(missing_ok=True)

        self.compressed = True

    def _get_h5_handle(self):
        if self.compressed:
            return h5py.File(self.compressed_file, "r")
        else:
            raise ValueError("Not compressed")

    def _get_full_data_dict(self, seed=None, sample_num=None):
        if self.compressed:
            raise ValueError

        seed_samplenum_str = self._seed_samplenum_str(seed, sample_num)

        if seed_samplenum_str == "":
            path = self._best_full_data_path

        else:
            path = self.dir_path / seed_samplenum_str / "confidences.json"

        return orjson.loads(open(path, "r").read())

    def _get_summary_data_dict(self, seed=None, sample_num=None):
        if self.compressed:
            raise ValueError

        seed_samplenum_str = self._seed_samplenum_str(seed, sample_num)

        if seed_samplenum_str == "":
            path = self._best_summary_path
        else:
            path = (
                self.dir_path / seed_samplenum_str / "summary_confidences.json"
            )

        return orjson.loads(open(path, "r").read())

    def _seed_samplenum_str(self, seed, sample_num):
        if self.server:
            raise NotImplementedError

        elif seed is None and sample_num is None:
            return ""

        sample_num = 0 if sample_num is None else sample_num

        if seed is None:
            raise ValueError("seed must be provided if sample_num is provided")

        return f"seed-{seed}_sample-{sample_num}"

    def _append_gz(self, path):
        return path.parent / (path.name + ".gz")
