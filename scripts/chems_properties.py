import os

from chems_pubchem_parse import ChemsParsePubchem


class ChemsProperties(ChemsParsePubchem):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_properties_dir = os.path.join(self.data_dir, 'chems_properties')
        self.chems_properties_assets_dir = os.path.join(self.data_dir, 'assets', 'chems_properties')