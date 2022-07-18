import numpy as np
from openbabel import pybel


def get_grid_center(ligand_pdbqt_filename: str, weighted_center: bool = False) -> np.ndarray:

    ligand = next(pybel.readfile(filename=ligand_pdbqt_filename, format="pdbqt"))
    coords = []
    for atom in ligand.atoms:
        coords.append(atom.coords)
    coords_array = np.array(coords)

    if weighted_center:
        center = np.mean(coords_array, axis=0)

    else:
        max_coords = np.amax(coords_array, axis=0)
        min_coords = np.amin(coords_array, axis=0)
        center = (max_coords + min_coords) / 2

    return center


def get_molecule_xyz(ligand_pdbqt_filename: str) -> np.ndarray:

    ligand = next(pybel.readfile(filename=ligand_pdbqt_filename, format="pdbqt"))
    coords = []
    for atom in ligand.atoms:
        coords.append(atom.coords)
    coords_array = np.array(coords)

    return coords_array


def cal_rmsd(molecule_1_coords: np.ndarray, molecule_2_coords: np.ndarray) -> float:

    assert molecule_1_coords.shape == molecule_2_coords.shape
    diff_sq = np.square(molecule_1_coords - molecule_2_coords)
    xyz_diff_sum = np.sum(diff_sq, axis=1)
    rmsd = np.sqrt(np.mean(xyz_diff_sum))

    return rmsd


def pdbqt_to_sdf(pdbqt_filename: str, output_filename: str) -> None:

    results = [m for m in pybel.readfile(filename=pdbqt_filename, format="pdbqt")]
    out = pybel.Outputfile(filename=output_filename, format="sdf", overwrite=True)
    for pose in results:

        pose.data.update({"Pose": pose.data["MODEL"]})
        pose.data.update({"Score": pose.data["REMARK"].split()[2]})
        del pose.data["MODEL"], pose.data["REMARK"], pose.data["TORSDO"]

        out.write(pose)
    out.close()
