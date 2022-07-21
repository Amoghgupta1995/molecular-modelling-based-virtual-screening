from openbabel import pybel


def sanitize_and_convert_ligand(
    input_ligand_filename: str, ligand_format: str, output_ligand_pdbqt_filename: str, overwrite: bool = True
) -> pybel.Molecule:
    ligand = next(pybel.readfile(filename=input_ligand_filename, format=ligand_format))
    ligand.addh()
    out = pybel.Outputfile(filename=output_ligand_pdbqt_filename, format="pdbqt", overwrite=True)
    out.write(ligand)
    out.close()

    return ligand
