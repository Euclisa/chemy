import os
from rdkit import Chem

from chems_reaction_properties import ChemsReactionProperties


class ChemsThermo(ChemsReactionProperties):
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.thermo_dir = os.path.join(self.data_dir, 'thermo')
    

    def _compute_formation_value(self, cid, value, cid_to_value):
        chem = self.cid_chem_map.get(cid)
        if chem is None:
            return None

        mol = Chem.MolFromInchi(chem['inchi'])
        if not mol:
            return None

        mol = Chem.AddHs(mol)
        
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        formation_value = value
        for symbol, count in atom_counts.items():
            el_entry = self.symb_to_el[symbol]
            el_cid = el_entry['cid']
            el_atom_count = el_entry['atom_count']
            el_thermo = cid_to_value.get(el_cid, None)

            if el_thermo is None:
                return None

            formation_value -= el_thermo * count / el_atom_count
        
        return formation_value