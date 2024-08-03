from device import setup_device
import numpy as np
from ase import io
from ase.build import fcc111
import math as m


def update_cell(system, cLength):
    a = system.cell[0]
    b = system.cell[1]
    c = np.array([0, 0, cLength])
    new_cell = np.array([a, b, c])
    system.set_cell(new_cell)
    mandev = setup_device(atoms=system)
    return mandev.shift_to_bottom()


# left electrode
spin = 'is'
ori_device = io.read('device-%s.vasp' % spin)
a_graphene = 2.467  # found by LDA relaxation
a_Ir = a_graphene*9./8.*m.sqrt(2)
Ir111 = fcc111('Ir', size=(1, 1, 3), a=a_Ir, vacuum=1.0, periodic=True)
mandev = setup_device(atoms=Ir111)
c_Ir = (mandev.get_zmax() - mandev.get_zmin()) * 3. / 2.
prin_layer = update_cell(Ir111, c_Ir)
io.write('left-electr.vasp', prin_layer)

# right electrode
Au_tip = ori_device[-12:].copy()
c_Au = (ori_device[-1].position[2]-ori_device[-4].position[2])*4
au_pl = update_cell(Au_tip, c_Au)
io.write('right-electr.vasp', au_pl)
