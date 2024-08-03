# build the scattering region: Autip--fepc--Gr9x9/Ir8x8
import numpy as np
from ase import io
from ase.build import fcc111
import math as m
from device import setup_device, setup_gr

a_graphene = 2.467  # found by LDA relaxation
d_Ir_gr = 3.23        # opt by LDA
d_gr_fepc = 3.20      # opt by LDA
d_fepc_tip = 5.0      # guess
a_Ir = a_graphene*9./8.*m.sqrt(2)

Ir111 = fcc111('Ir', size=(8, 8, 6), a=a_Ir, vacuum=17.5, periodic=True)
del Ir111[-8*8*2:]
dev = setup_device(atoms=Ir111)
dev.shift_to_bottom()
gr = setup_gr(a=a_graphene, gamma=60, size=(9, 9, 1))
gr.translate([-gr.cell[1]])
dev.add_atoms(in_atoms=gr, offset=d_Ir_gr, atoms_pbc=True)
tip = io.read('tip.xyz')
for spin in ['hs', 'is']:
    device = dev.copy()
    molecule = io.read('fe-%s-LDA-U4.0.STRUCT_OUT' % spin)
    device.add_atoms(in_atoms=molecule, offset=d_gr_fepc, atoms_pbc=False)
    tip.set_cell(device.atoms.cell)
    tip.center()
    aupos = tip.get_positions()
    tip_atom_index = np.argmin(aupos[:, 2])
    pos_Au = aupos[tip_atom_index]
    for atom in device.atoms:
        if atom.symbol == 'Fe':
            pos_Fe = atom.position
    vector = pos_Fe-pos_Au
    shift = [vector[0], vector[1], 0]
    tip.translate(shift)        # shift the tip to point at Fe
    device.add_atoms(in_atoms=tip, offset=d_fepc_tip, atoms_pbc=True)
    io.write('device-%s.vasp' % spin, device.atoms)
