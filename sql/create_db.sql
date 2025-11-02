SET search_path TO public;


CREATE TABLE compounds (
    cid INTEGER PRIMARY KEY,
    name VARCHAR(1000) NOT NULL,
    mf VARCHAR(100) NOT NULL,
    mw REAL NOT NULL,
    charge SMALLINT NOT NULL,
    smiles VARCHAR(1000) NOT NULL,
    inchi VARCHAR(1000) NOT NULL,
    inchikey CHAR(27) UNIQUE NOT NULL,
    complexity REAL,
    bertz_complexity REAL NOT NULL,
    organic BOOLEAN NOT NULL
);


CREATE TABLE compound_cas (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    cas VARCHAR(15) NOT NULL,
    PRIMARY KEY (cid, cas)
);


CREATE TABLE compound_fingerprints (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    ECFP4_fp INTEGER[64] NOT NULL,
    popcount SMALLINT NOT NULL
);


CREATE TABLE compound_synonyms (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    synonym VARCHAR(1000) NOT NULL,
    PRIMARY KEY (cid, synonym)
);


CREATE TABLE compound_wiki (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    wiki VARCHAR(300) NOT NULL
);


CREATE TABLE compound_nfpa (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    health REAL,
    flammability REAL,
    instability REAL
);


CREATE TABLE categories (
    code VARCHAR(10) PRIMARY KEY,
    name VARCHAR(100) NOT NULL,
    domain VARCHAR(10) NOT NULL
);


CREATE TABLE compound_categories (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    category_code VARCHAR(10) NOT NULL REFERENCES categories(code),
    PRIMARY KEY (cid, category_code)
);


CREATE TABLE compound_hazard_pictograms (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    pictogram CHAR(5) NOT NULL,
    PRIMARY KEY (cid, pictogram)
);


CREATE TABLE compound_properties (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    property_name VARCHAR(100) NOT NULL,
    property_value VARCHAR(3000) NOT NULL,
    PRIMARY KEY (cid, property_name)
);


CREATE TABLE compound_commonness_sorting (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    rank INTEGER UNIQUE NOT NULL
);


CREATE TABLE compound_descriptions (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    description TEXT NOT NULL
);


CREATE TABLE compound_complexity_sorting (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    rank INTEGER UNIQUE NOT NULL
);

CREATE TABLE compound_curiosity_sorting (
    cid INTEGER PRIMARY KEY REFERENCES compounds(cid),
    rank INTEGER UNIQUE NOT NULL
);


CREATE TABLE reactions (
    rid CHAR(24) PRIMARY KEY,
    complexity REAL NOT NULL,
    source VARCHAR(100) NOT NULL,
    balanced BOOLEAN NOT NULL,
    confidence REAL
);


CREATE TABLE reaction_reactants (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    coeff INTEGER CHECK (coeff > 0 OR coeff IS NULL),
    rid CHAR(24) NOT NULL REFERENCES reactions(rid),
    PRIMARY KEY (cid, rid)
);


CREATE TABLE reaction_products (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    coeff INTEGER CHECK (coeff > 0 OR coeff IS NULL),
    rid CHAR(24) NOT NULL REFERENCES reactions(rid),
    PRIMARY KEY (cid, rid)
);


CREATE TABLE reaction_solvents (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    rid CHAR(24) NOT NULL REFERENCES reactions(rid),
    PRIMARY KEY (cid, rid)
);


CREATE TABLE reaction_catalysts (
    cid INTEGER NOT NULL REFERENCES compounds(cid),
    rid CHAR(24) NOT NULL REFERENCES reactions(rid),
    PRIMARY KEY (cid, rid)
);


CREATE TABLE reaction_details (
    rid CHAR(24) PRIMARY KEY REFERENCES reactions(rid),
    doi VARCHAR(100),
    patent VARCHAR(100),
    description TEXT,
    source VARCHAR(100) NOT NULL,
    confidence REAL
);