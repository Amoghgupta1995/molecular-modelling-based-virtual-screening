from openmm.app import PDBFile
from pdbfixer import PDBFixer


def fix_pdb(
    input_pdb_filename: str,
    output_pdb_filename: str,
    find_missing_residues: bool = True,
    find_replace_non_standard_residues: bool = True,
    keep_waters: bool = True,
    find_add_missing_heavy_atoms: bool = True,
    protonation_pH: float = 7.0,
) -> None:

    fixer = PDBFixer(filename=input_pdb_filename)

    if find_missing_residues:
        fixer.findMissingResidues()

    if find_replace_non_standard_residues:
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

    if keep_waters:
        fixer.removeHeterogens(True)
    else:
        fixer.removeHeterogens(False)

    if find_add_missing_heavy_atoms:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

    fixer.addMissingHydrogens(protonation_pH)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_filename, "w"))
