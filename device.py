from ase.lattice.hexagonal import Graphene
import numpy as np


def atoms_separation(bulk):
    pos = bulk.get_positions()
    zm = max(pos[:, 2])
    shift = bulk.cell[2, 2] - zm
    return shift


class setup_device:
    """Attach electrodes to device region"""

    def __init__(self, atoms=None):
        self.atoms = atoms.copy()

    def get_zmin(self):
        pos = self.atoms.get_positions()
        return min(pos[:, 2])

    def get_zmax(self):
        pos = self.atoms.get_positions()
        return max(pos[:, 2])

    def move(self, s=0.0):
        for atom in self.atoms:
            atom.z += s

    def expand_cell(self, in_atoms=None):
        a = self.atoms.cell[0]
        b = self.atoms.cell[1]
        c = self.atoms.cell[2] + in_atoms.cell[2]
        new_cell = np.array([a, b, c])
        self.atoms.set_cell(new_cell)

    def shift_c(self, c=0.0):
        self.move(c)
        return self.atoms

    def shift_to_bottom(self):
        zmin = self.get_zmin()
        self.move(s=-zmin)
        return self.atoms

    def add_atoms(self, in_atoms=None, offset=0.0, atoms_pbc=False):
        """The crystal cell is fixed.
        When adding periodic structures (atoms_pbc=True)
        the lattice in xy-plane must be identical to self.atoms"""
        if not atoms_pbc:
            in_atoms.set_cell(self.atoms.cell)
            in_atoms.center()
        zmax = self.get_zmax()
        shift = zmax + offset
        layer = self.__class__(atoms=in_atoms)
        layer.shift_to_bottom()
        layer.shift_c(c=shift)
        self.atoms.extend(layer.atoms)
        return self.atoms

    def copy(self):
        return self.__class__(self.atoms)

    def attach_electrodes(self,
                          left_electr=None,
                          right_electr=None):
        """Cell must match in xy-plane"""
        # attach left electrode
        if left_electr is not None:
            sl = atoms_separation(bulk=left_electr)
            electrL = self.__class__(atoms=left_electr)
            electrL.add_atoms(in_atoms=self.atoms, offset=sl, atoms_pbc=True)
            electrL.expand_cell(in_atoms=self.atoms)
            self.atoms = electrL.atoms

        # attach right electrode
        if right_electr is not None:
            sr = atoms_separation(bulk=right_electr)
            self.add_atoms(in_atoms=right_electr, offset=sr, atoms_pbc=True)
            self.expand_cell(in_atoms=right_electr)

        return self.atoms

    def get_geometry(self):
        return self.atoms


def setup_gr(a, gamma=120, size=(1, 1, 1), c=3.355, vacuum=0.0):
    gr = Graphene(symbol='C',
                  latticeconstant={
                      'a': a,
                      'c': c,
                  },
                  size=size,
                  pbc=(1, 1, 1))
    if abs(vacuum) > 1.e-3:
        gr.center(vacuum=vacuum, axis=2)
    if gamma < 90:
        a = gr.cell[0]
        b = np.array([gr.cell[1, 0], -gr.cell[1, 1], gr.cell[1, 2]])
        c = gr.cell[2]

        new_cell = np.array([a, b, c])
        gr.set_cell(new_cell, scale_atoms=True)
    return gr


def getMoleculeIndex(atoms, gr=(9, 9), chem=['C', 'H', 'N', 'Fe'],
                     hasGraphene=True):
    index = []
    for atom in atoms:
        for element in chem:
            if atom.symbol == element:
                index.append(atom.index)
                break
    if hasGraphene:
        if 'C' in chem:
            nC_in_gr = gr[0]*gr[1]*2
            del index[:nC_in_gr]
    return index


def update_cell(system, cLength):
    a = system.cell[0]
    b = system.cell[1]
    c = np.array([0, 0, cLength])
    new_cell = np.array([a, b, c])
    system.set_cell(new_cell)
    mandev = setup_device(atoms=system)
    return mandev.shift_to_bottom()


def get_local_spin(atoms, mag_centers=['fe'], spin_states=['ls'], s=None):
    spin = {
        'Mn': (5.0, 1.0),
        'Fe': (4.0, 0.0, 2.0),
        'Co': (3.0, 1.0),
        'Ni': (2.0, 2.0)
    }
    number_of_atoms = atoms.get_global_number_of_atoms()
    local_magetic_moment = np.zeros(number_of_atoms)
    for i, metal in enumerate(mag_centers):
        magcenter = metal.capitalize()
        if s is not None:
            local_spin = s[i]
        else:
            spinstate = spin_states[i].lower()
            if spinstate == 'hs':
                local_spin = spin[magcenter][0]
            elif spinstate == 'is':
                local_spin = spin[magcenter][2]
            else:
                local_spin = spin[magcenter][1]
        for atom in atoms:
            if atom.symbol == magcenter:
                local_magetic_moment[atom.index] = local_spin
    return local_magetic_moment
