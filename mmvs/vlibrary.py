from pathlib import Path
from openbabel import pybel


class VirtualLibrary:
    def __init__(self, virtual_library_name: str, virtual_library_path: Path):
        self.virtual_library_name = virtual_library_name
        if virtual_library_path.exists():
            self.virtual_library_path = virtual_library_path
            self.virtual_library_type = self.virtual_library_path.suffix
        else:
            raise FileNotFoundError("Virtual library file does not exist.")
        self.virtual_library = None

    def load(self):
        if self.virtual_library_path is not None:
            self.virtual_library = pybel.readfile(
                # Remove "." from ".sdf" and convert Path type to str as Openbabel expects string.
                format=self.virtual_library_type.replace(".", ""), filename=str(self.virtual_library_path)
            )
            for mol in self.virtual_library:
                print(type(mol))
