import numpy as np
import orjson
from pathlib import Path
import polars as pl
import MDAnalysis as mda
import h5py
import gzip

class BoltzOutput():

    def __init__(self, directory_path, **kwargs):

        dir_path = None
        if isinstance(directory_path, str):
            dir_path = Path(directory_path)
        elif isinstance(directory_path, Path):
            dir_path = directory_path
        else:
            raise ValueError

        self.job_name = dir_path.name
        self.dir_path = dir_path

        self.compressed_file = Path(self.dir_path / f"{self.job_name}.h5")

        if self.compressed_file.is_file():
            self.compressed = True
        else:
            self.compressed = False

        best_model = dir_path / f"{self.job_name}_model_0.cif"
        fmt = "MMCIF"
        if not best_model.exists():
            best_model = dir_path / f"{self.job_name}_model_0.pdb"
            fmt = "pdb"
            if not best_model.exists():
                raise FileExistsError(f"No structure file found for {directory_path}")

        self._fmt = fmt
    

    def has_affinity(self):
        return (self.dir_path/ f"affinity_{self.job_name}.json").exists()
    
    def get_affinity_metrics(self):
        raise NotImplementedError
    
    def get_summary_metrics(self, sample_num=0):
        
        if self.compressed:
            raise NotImplementedError
        
        else:
            summ_dict = orjson.loads(open(self.dir_path / f"confidence_{self.job_name}_model_{sample_num}.json"))
            return summ_dict
    
    def get_pae_ndarr(self, sample_num=0):
        
        if self.compressed:
            raise NotImplementedError
        
        return np.load(self.dir_path / f"pae_{self.job_name}_model_{sample_num}.npz")["pae"]

    def get_pde_ndarr(self, sample_num=0):
        if self.compressed:
            raise NotImplementedError

        return np.load(self.dir_path / f"pde_{self.job_name}_model_{sample_num}.npz")["pde"]
        
    def get_plddt_ndarr(self, sample_num=0):
        if self.compressed:
            raise NotImplementedError
        
        return np.load(self.dir_path / f"plddt_{self.job_name}_model_{sample_num}.npz")["plddt"]
        
    def get_mda_universe(self, sample_num=0, **kwargs):
        if self.compressed:
            raise NotImplementedError

        top_path = (self.dir_path / f"{self.job_name}_model_{sample_num}.{self.fmt}")


        return mda.Universe(top_path, topology_format=self.fmt, **kwargs)
        

    def compress(self):
        raise NotImplementedError    