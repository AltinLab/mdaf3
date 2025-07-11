from MDAnalysis.analysis.dssp import DSSP
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import polars as pl
from itertools import repeat


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
            curr_helix = [int(helix_resindices[i])]

        prev_ix = helix_resindices[i]

    return helices


def _apply_row(row, func, args, kwargs):
    return func(row, *args, **kwargs)


def split_apply_combine(df, func, *args, **kwargs):
    """
    Applies `func(row, *args, **kwargs)` to each row of `df` in parallel,
    then concatenates the returned DataFrames.
    """
    # an iterator over your rows
    rows = df.iter_rows(named=True)

    if "chunksize" in kwargs:
        chunksize = kwargs.pop("chunksize")
    else:
        chunksize = 1

    partials = process_map(
        _apply_row,
        rows,
        repeat(func),
        repeat(args),
        repeat(kwargs),
        chunksize=chunksize,
        total=df.height,
    )

    return pl.DataFrame(partials, infer_schema_length=len(partials))


def serial_apply(df, func, *args, **kwargs):

    out = []
    for row in tqdm(
        df.iter_rows(named=True), total=df.height, desc="Processing rows"
    ):
        out.append(func(row, *args, **kwargs))

    return pl.DataFrame(out, infer_schema_length=len(out))
