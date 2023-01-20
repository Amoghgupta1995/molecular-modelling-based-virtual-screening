import hashlib
from pathlib import Path

from loguru import logger
from openbabel import pybel
from tqdm import tqdm


class VirtualLibrary:
    def __init__(self, virtual_library_name: str, virtual_library_path: Path, output_directory: str):
        self.virtual_library_name = virtual_library_name
        self.current_working_directory = Path.cwd()
        if virtual_library_path.exists():
            self.virtual_library_path = virtual_library_path
            self.virtual_library_type = self.virtual_library_path.suffix
        else:
            raise FileNotFoundError("Virtual library file does not exist.")
        self.output_directory = output_directory
        self.virtual_library = None
        self.virtual_library_df = None
        self.mol_pdbqt_dict: dict[str, str] = dict()
        logger.success("Virtual library successfully initialized.")

    def load(self):

        pybel.ob.obErrorLog.SetOutputLevel(0)

        smiles_list = []
        smiles_md5_hash_list = []
        molecular_weight = []

        if self.virtual_library_path is not None:
            self.virtual_library = list(pybel.readfile(
                # Remove "." from ".sdf" and convert Path type to str as Openbabel expects string.
                format=self.virtual_library_type.replace(".", ""),
                filename=str(self.virtual_library_path),
            ))

        logger.info("Parsing, encoding and converting the format of molecules in the Virtual library.")

        molecule: pybel.Molecule
        for molecule in tqdm(self.virtual_library, colour="green"):
            smiles = str(molecule.write(format="smi"))
            smiles_md5_hash = hashlib.md5(smiles.encode("utf-8")).hexdigest()
            smiles_list.append(smiles)
            smiles_md5_hash_list.append(smiles_md5_hash)
            molecule.make3D()
            molecular_weight.append(molecule.molwt)
            self.mol_pdbqt_dict[smiles_md5_hash] = molecule.write(format="pdbqt")

        test_pdbqt = list(self.mol_pdbqt_dict.values())[10]
        with open("test.pdbqt", 'w') as f:
            f.write(test_pdbqt)
