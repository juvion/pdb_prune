"""
Convert NumPy coordinate arrays back to minimal PDB format.
"""

import numpy as np
from pathlib import Path


class NumpyToPDBConverter:
    """Re-create a readable PDB file from a NumPy coordinate array."""

    def __init__(self, atom_names=None, res_names=None, res_ids=None, chains=None):
        """
        Optionally supply fixed atom / residue metadata.
        If not provided, defaults are generated.
        """
        self.atom_names = atom_names or ["C", "N", "O", "P"]
        self.res_names  = res_names  or ["A", "U", "G", "C"]
        self.res_ids    = res_ids
        self.chains     = chains

    # ------------------------------------------------------------------ #
    # Public API                                                         #
    # ------------------------------------------------------------------ #

    def convert(self, coords, out_pdb, res_ids=None, chains=None, atom_names=None):
        """
        coords : (N, 3) NumPy array of XYZ coordinates
        out_pdb: pathlib.Path or str – file to write
        res_ids: iterable length N (1-based residue numbers)
        chains : iterable length N (chain IDs)
        atom_names: iterable length N (atom names)
        """
        coords   = np.asarray(coords, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError("coords must be (N, 3)")

        n_atoms = coords.shape[0]

        # Defaults if caller didn’t supply per-atom metadata
        res_ids    = res_ids    or (np.arange(n_atoms) // 4 + 1)
        chains     = chains     or (["A"] * n_atoms)
        atom_names = atom_names or (["P", "C1'", "N1", "C2"] * (n_atoms // 4 + 1))[:n_atoms]

        lines = []
        for i, (x, y, z) in enumerate(coords, start=1):
            lines.append(
                f"ATOM  {i:5d} {atom_names[i-1]:<4} {self.res_names[0]:3} "
                f"{chains[i-1]:1}{res_ids[i-1]:>4}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        lines.append("END\n")

        Path(out_pdb).write_text("".join(lines))
        return out_pdb