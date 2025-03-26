mdaf3
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  |
| :----------------- | :------- |
| **Status**         | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Community**      | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

[badge_actions]: https://github.com/ljwoods2/mdaf3/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/ljwoods2/mdaf3/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/ljwoods2/mdaf3/latest
[badge_docs]: https://readthedocs.org/projects/mdaf3/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/ljwoods2/mdaf3.svg
[url_actions]: https://github.com/ljwoods2/mdaf3/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/ljwoods2/mdaf3/branch/main
[url_docs]: https://mdaf3.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/ljwoods2/mdaf3/releases
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org

> [!IMPORTANT]
> This package relies on a feature branch of MDAnalysis to parse mmCIF files and
> therefore should be considered experimental. In the future, if this feature is merged,
> this package will use the official release. For now, please carefully validate your
> selections and see [tests](mdaf3/tests/test_mdaf3.py) to see what is currently validated.

[AlphaFold3](https://github.com/google-deepmind/alphafold3) outputs a set of [confidence metrics](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md) that are useful for i.e. protein binding predictions, however, making use of these metrics requires careful and time consuming parsing.

[MDAnalysis](https://www.mdanalysis.org/) provides an [atom selection language](https://userguide.mdanalysis.org/stable/selections.html) (think SQL for molecular topologies) that makes associating confidence metrics with molecular positions, amino acid/atom type, and other topological information easy.

This package seeks to expose AF3 outputs (including confidence metrics) and predicted topologies via an easy-to-use interface.

### How do I...

Get a python interface into my AF3 output?
```python
# example inference output
from mdaf3.data.files import UNCOMPRESSED_AF3_OUTPUT_PATH
from mdaf3.AF3OutputParser import AF3Output

# Equivalent to AF3Output("/path/to/inference/output")
af3_output = AF3Output(UNCOMPRESSED_AF3_OUTPUT_PATH)
```

Get summary information about the AF3 inference job?
```python
summary_dict = af3_output.get_summary_metrics()

chain_pair_iptm_ndarr = summary_dict["chain_pair_iptm"]
ranking_score = summary_dict["ranking_score"]

# all 'get_' methods take an optional "seed" and "sample_num" argument
# if not provided, the best model (by AF3 ranking score)
# is returned
summary_dict_seed_1 = af3_output.get_summary_metrics(seed=1)
```

Calculate the mean pLDDT of atoms in a particular protein?
```python
u = af3_output.get_mda_universe()

# segid in MDAnalysis corresponds to protein ID in AF3 input JSON
# https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
selection = u.select_atoms("segid A")

# pLDDT is stored in the "tempfactors" attribute of an 
# MDAnalysis AtomGroup
mean_pLDDT_chain = selection.tempfactors.mean()
```

Find the minimum PAE among all residue pairs between two protein chains?
```python
u = af3_output.get_mda_universe()

# move up the topology heirarchy
# from atom -> amino acid residue
# using AtomGroup.residues
protein_1_all_res = u.select_atoms("segid A").residues

protein_2_all_res = u.select_atoms("segid B").residues

pae_ndarr = af3_output.get_pae_ndarr()

# resindices are 0-indexed amino acid residue indices
# that correspond to AF3 token indices
min_pae_p1_p2 = pae_ndarr[protein_1_res.resindices][:,
    protein_2_all_res.resindices
].min()
```

Find the max contact probability between a single residue of one protein chain
and any residue in another protein chain?
```python
u = af3_output.get_mda_universe()

protein_res_1 = u.select_atoms("segid A").residues[0]
protein_2_all_res = u.select_atoms("segid B").residues

contact_prob_ndarr = af3_output.get_contact_prob_ndarr()

max_contact_prob_res1_p2 = contact_prob_ndarr[protein_res_1.resindices][:,
    protein_2_all_res.resindices
].max()
```

Find the mean pLDDT of all atoms that are within 5 angstroms 
of a particular residue?
```python
u = af3_output.get_mda_universe(seed=1)

particular_residue = u.select_atoms("segid A").residues[0]

atoms_around_particular_res = u.select_atoms(
    "around 5 group pr", pr=particular_residue
)

mean_pLDDT_around_pr = atoms_around_particular_res.mean()
```

Batch apply a feature extraction method to all my AF3 jobs (with job names
stored in a Polars DataFrame)
```python
from pathlib import Path
import polars as pl
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import serial_apply, split_apply_combine

def extract_protein1_mean_pLDDT(row, af3_parent_dir):
    job_dir = Path(af3_parent_dir) / row["job_name"]
    af3_output = AF3Output(job_dir)
    u = af3_output.get_mda_universe()
    protein1_mean_pLDDT = u.select_atoms("segid A").tempfactors.mean()
    row["protein1_mean_pLDDT"] = protein1_mean_pLDDT
    return pl.DataFrame(row)

all_jobs = pl.DataFrame({"job_name": ["93f0240a1d2c15da9551841d22239d41"]})

af3_parent_dir = "mdaf3/data"

# use split_apply_combine for process-parallel execution.
# these methods will convert each row into a dict,
# pass it to your extraction method,
# and then concat the resulting pl.DataFrames
all_job_with_feat = serial_apply(
    all_jobs, extract_protein1_mean_pLDDT, af3_parent_dir
)

feature_np = (
    all_job_with_feat.select("protein1_mean_pLDDT").to_series().to_numpy()
)
```

Compress my AF3 output directory without losing confidence metric precision?

> [!NOTE]
> This will delete 'TERMS_OF_USE.md' as well as the input JSON for the AF3 job ('<job_name>_data.json')
> among other things. This feature is designed with large HPC batches in mind, so if you aren't sure, 
> read the [compression code](https://github.com/ljwoods2/mdaf3/blob/main/mdaf3/AF3OutputParser.py#:~:text=compress)!

```python
af3_output = AF3Output(UNCOMPRESSED_AF3_OUTPUT_PATH)
af3_output.compress()
```


### Installation

Below we provide instructions both for `conda` and
for `pip`.

First, clone the repo locally:

```
git clone https://github.com/ljwoods2/mdaf3.git
cd mdaf3
```

#### With conda

Ensure that you have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

Create a virtual environment and activate it:

```
conda create --name mdaf3
conda activate mdaf3
```

Install the dependencies:

```
conda env update --name mdaf3 --file devtools/conda-envs/test_env.yaml
```

Build this package from source:

```
pip install -e .
```

If you want to update your dependencies (which can be risky!), run:

```
conda update --all
```

And when you are finished, you can exit the virtual environment with:

```
conda deactivate
```

#### With pip

To build the package from source, run:

```
pip install .
```

If you want to create a development environment, install
the dependencies required for tests and docs with:

```
pip install ".[test]"
```

mdaf3 is bound by a [Code of Conduct](https://github.com/ljwoods2/mdaf3/blob/main/CODE_OF_CONDUCT.md).

### Copyright

The mdaf3 source code is hosted at https://github.com/ljwoods2/mdaf3
and is available under the GNU General Public License, version 2 (see the file [LICENSE](https://github.com/ljwoods2/mdaf3/blob/main/LICENSE)).

Copyright (c) 2025, Lawson Woods


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using mdaf3 in published work.
