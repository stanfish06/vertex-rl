"""Build and run the focused native mesh and contact checks in unittest discovery."""
from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "model" / "three_d"


class Native3D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        cls.scratch = Path(workspace.name)
        cls.compiler = [os.environ.get("CXX", "g++"), "-std=c++17", "-O2", "-I", str(MODEL)]
        cls.mesh_object = cls.scratch / "mesh.o"
        subprocess.run([*cls.compiler, "-c", MODEL / "mesh.cpp", "-o", cls.mesh_object], check=True)

    def test_mesh_geometry(self):
        binary = self.scratch / "test-mesh"
        subprocess.run([*self.compiler, ROOT / "tests" / "mesh3d.cpp",
                        self.mesh_object, "-o", binary], check=True)
        subprocess.run([binary, self.scratch / "imports"], check=True)

    def test_contact_geometry(self):
        binary = self.scratch / "test-contact"
        subprocess.run([*self.compiler, ROOT / "tests" / "test_contact3d.cpp",
                        MODEL / "contact.cpp", self.mesh_object, "-o", binary], check=True)
        subprocess.run([binary], check=True)


if __name__ == "__main__":
    unittest.main()
