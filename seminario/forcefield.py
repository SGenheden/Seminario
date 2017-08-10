"""
This module contains tools to estimate force field parameters with the
Seminario method.
"""
from __future__ import division, print_function, absolute_import

import numpy as np
import parmed

__all__ = ['seminario_bond', 'make_bonded_ff']

def seminario_bond(bond, hessian, scaling) :
    """
    Estimate the bond force constant using the Seminario method, i.e. by
    analysing the Hessian submatrix. Will average over atom1->atom2 and
    atom2->atom1 force constants

    Parameters
    ----------
    bond : parmed.Bond
        the bond to estimate the force constant for
    hessian : numpy.ndarray
        the Hessian matrix
    scaling : float
        the Hessian scaling factor
    """

    def _calc_force(atm1, atm2, hess, vec12):
        submat = -hess[3*atm1.idx:3*atm1.idx+3, 3*atm2.idx:3*atm2.idx+3]
        eigval, eigvec = np.linalg.eig(submat)
        fc = 0.0
        for i in range(3):
            fc += eigval[i] * np.abs(np.dot(eigvec[:,i], vec12))
        return fc

    atm1 = bond.atom1
    atm2 = bond.atom2

    # Vector from atm1 to atm2
    vec12 = np.asarray([atm1.xx - atm2.xx, atm1.xy - atm2.xy, atm1.xz - atm2.xz])
    vec12 /= bond.measure()

    f12 = _calc_force(atm1, atm2, hessian, vec12)
    f21 = _calc_force(atm2, atm1, hessian, vec12)
    # 2240.87 is from Hartree/Bohr ^2 to kcal/mol/A^2
    # 418.4 is kcal/mol/A^2 to kJ/mol/nm^2
    return scaling * 2240.87 * 418.4 * 0.5 * (f12+f21)

def make_bonded_ff(struct, xyz_orig, hess, bonds, scaling) :
    """
    Make bonded force field parameters for selected bonds using the Seminario
    method.

    This will print the force field parameters to standard output

    Parameters
    ----------
    struct : parmed.Structure
        the structure of the model, including optimized coordinates
    xyz_orig : numpy.ndarray
        the original xyz coordinates
    hess : numpy.ndarray
        the Hessian matrix
    bonds : list of strings
        the bond definitions
    scaling : float
        the scaling factor for the Hessian
    """

    for bond in bonds :
        sel1, sel2 = bond.strip().split("-")
        atom1 = struct.view[sel1].atoms[0]
        atom2 = struct.view[sel2].atoms[0]
        bond = parmed.Bond(atom1, atom2)
        struct.bonds.append(bond)
    struct_orig = struct.copy(type(struct))
    struct_orig.coordinates = xyz_orig

    print("Bond\tidx1\tidx2\tr(x-ray)\tr(opt) [nm]\tk [kJ/mol/nm2]")
    for bond, bond_orig, bondmask in zip(struct.bonds, struct_orig.bonds, bonds) :
        force = seminario_bond(bond, hess, scaling)
        print("%s\t%d\t%d\t%.4f\t%.4f\t%.4f"%(bondmask, bond.atom1.idx+1,
            bond.atom2.idx+1, bond_orig.measure()*0.1, bond.measure()*0.1, force))
