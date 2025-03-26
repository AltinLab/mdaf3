from MDAnalysis.analysis.dssp import DSSP
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import polars as pl


def raw_helix_indices(sel):
    # find helices
    # https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/dssp.html
    helix_resindices_boolmask = DSSP(sel).run().results.dssp_ndarray[0, :, 1]
    return sel.residues[helix_resindices_boolmask].resindices


def anneal_helix_indices(helix_resindices, len_tol=15, gap_tol=10):
    """
    Parameters
    ----------
    helix_resindices : list
        List of residue indices that are part of a helix.
    len_tol : int, optional
        Minimum length of a helix. The default is 15.
    gap_tol : int, optional
        Maximum gap between residues in a helix. The default is 10.
    """

    helices = []
    curr_helix = []
    prev_ix = helix_resindices[0] - 1
    for i in range(len(helix_resindices)):

        if helix_resindices[i] - prev_ix <= gap_tol:
            # fill gap
            append_ix = range(prev_ix + 1, helix_resindices[i] + 1)
            curr_helix.extend(append_ix)

            if i == len(helix_resindices) - 1 and len(curr_helix) >= len_tol:
                helices.append(curr_helix)

        else:
            if len(curr_helix) >= len_tol:
                helices.append(curr_helix)
            # start new helix
            curr_helix = [helix_resindices[i]]

        prev_ix = helix_resindices[i]

    return helices


def split_apply_combine(df, func, *args, **kwargs):
    results = []
    # Create a list of rows to process.
    rows = list(df.iter_rows(named=True))

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(func, row, *args, **kwargs) for row in rows]
        for future in as_completed(futures):
            results.append(future.result())

    return pl.concat(results)


def serial_apply(df, func, *args, **kwargs):

    out = []
    for row in tqdm(
        df.iter_rows(named=True), total=df.height, desc="Processing rows"
    ):
        out.append(func(row, *args, **kwargs))

    return pl.concat(out)
