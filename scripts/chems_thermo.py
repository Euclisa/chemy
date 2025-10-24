import os

from chems_parse_reactions import ChemsParseReactions


class ChemsThermo(ChemsParseReactions):
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.thermo_dir = os.path.join(self.data_dir, 'thermo')