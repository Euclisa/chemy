from scripts.chemy import Chemy


def test_chemy_wires_composed_services(sample_data_dir):
    chemy = Chemy(str(sample_data_dir))

    assert chemy.compiler is chemy.properties
    assert chemy.solubility.crc is chemy.crc
    assert chemy.properties.reactions is chemy.reactions
    assert chemy.llm.llm_client is chemy.llm_client
    assert chemy.sql.properties is chemy.properties
    assert chemy.thermo_burcat.thermo is chemy.thermo
    assert chemy.thermo_experiments.properties is chemy.properties
    assert chemy.thermo_xtb.thermo is chemy.thermo
    assert chemy.thermo_llm.reaction_parser is chemy.reaction_llm
