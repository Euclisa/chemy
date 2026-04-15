import os

from rdkit import Chem


class ThermoOps:
    def __init__(self, data_dir, compounds):
        self.data_dir = data_dir
        self.compounds = compounds

        self.thermo_dir = os.path.join(self.data_dir, 'thermo')

    def get_atom_counts(self, cid):
        chem = self.compounds.cid_chem_map.get(cid)
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

        return atom_counts

    def get_element_reference_coeffs(self, cid):
        atom_counts = self.get_atom_counts(cid)
        if atom_counts is None:
            return None

        coeffs = {}
        for symbol, count in atom_counts.items():
            el_entry = self.compounds.symb_to_el.get(symbol)
            if el_entry is None:
                return None

            el_cid = el_entry['cid']
            el_atom_count = el_entry['atom_count']
            coeffs[el_cid] = coeffs.get(el_cid, 0) + count / el_atom_count

        return coeffs

    def compute_formation_value(self, cid, value, cid_to_value):
        element_coeffs = self.get_element_reference_coeffs(cid)
        if element_coeffs is None:
            return None

        formation_value = value
        for el_cid, coeff in element_coeffs.items():
            el_thermo = cid_to_value.get(el_cid)
            if el_thermo is None:
                return None

            formation_value -= el_thermo * coeff

        return formation_value

    def compute_reaction_value(self, reaction, cid_to_value, balance):
        result = 0.0
        for entry in reaction['reagents']:
            value = cid_to_value.get(entry['cid'])
            if value is None:
                return None
            result -= balance[entry['cid']] * value

        for entry in reaction['products']:
            value = cid_to_value.get(entry['cid'])
            if value is None:
                return None
            result += balance[entry['cid']] * value

        return result
