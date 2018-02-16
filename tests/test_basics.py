
from __future__ import division, print_function, absolute_import

import os
import unittest

import seminario
import parmed
import numpy as np

fchk_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..","examples","model.fchk")
gro_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..","examples","model.gro")


class TestSeminario(unittest.TestCase):

    def test_loadstruct(self) :

        struct = parmed.load_file(gro_filename, skip_bonds=True)
        self.assertEqual(len(struct.atoms),41)

    def test_loadfchk(self) :

        xyz_opt, hess = seminario.gaussian.parse_fchk(fchk_filename)
        self.assertEqual(xyz_opt.shape,(41,3))
        self.assertEqual(hess.shape,(123,123))
        self.assertEqual(xyz_opt[0,0],-2.39626372)
        self.assertEqual(hess[0,0],0.493623378)

    def test_bondparm(self) :

        struct = parmed.load_file(gro_filename, skip_bonds=True)
        xyz_opt, hess = seminario.gaussian.parse_fchk(fchk_filename)
        struct.coordinates = xyz_opt*0.529177249
        atom1 = struct.view[":1@N"].atoms[0]
        atom2 = struct.view[":3@CU"].atoms[0]
        bond = parmed.Bond(atom1, atom2)
        force = seminario.forcefield.seminario_bond(bond, hess, 0.963)
        self.assertEqual(np.round(force,4), 60493.9023)
        self.assertEqual(np.round(bond.measure(),3), 1.993)

    def test_angleparm(self) :

        struct = parmed.load_file(gro_filename, skip_bonds=True)
        xyz_opt, hess = seminario.gaussian.parse_fchk(fchk_filename)
        struct.coordinates = xyz_opt*0.529177249
        atom1 = struct.view[":1@CG"].atoms[0]
        atom2 = struct.view[":1@ND1"].atoms[0]
        atom3 = struct.view[":3@CU"].atoms[0]
        angle = parmed.Angle(atom1, atom2, atom3)
        force = seminario.forcefield.seminario_angle(angle, hess, 0.963)
        self.assertEqual(np.round(force,4), 192.3902)
        self.assertEqual(np.round(angle.measure(),3), 123.642)

if __name__ == '__main__':
    
    unittest.main()
