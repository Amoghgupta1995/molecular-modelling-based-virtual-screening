from typing import List, Tuple

from openbabel import pybel
from tqdm import tqdm


class VirtualScreen(object):
    def __init__(self, experiment_name: str):
        self.experiment_name = experiment_name

    def load_ligand_library(self, sdf_file_path: str):
        self.ligand_library = pybel.readfile(format="sdf", filename=sdf_file_path)

    def show_molecule_weights(self):
        if self.ligand_library is not None:
            for mol in self.ligand_library:
                mol.addh()
                if mol is not None:
                    print(mol.molwt)

    def show_library_size(self):
        print(f"Library size is {len(list(self.ligand_library))}.")

