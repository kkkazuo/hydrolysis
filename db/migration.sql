CREATE TABLE status (
    id VARCHAR(68) NOT NULL UNIQUE,
    status INTEGER NOT NULL
);

CREATE TABLE result (
    id VARCHAR(68) NOT NULL UNIQUE,
    chem_index INTEGER NOT NULL,
    smiles TEXT NOT NULL,
    hydrolyzability INTEGER,
    FOREIGN KEY(id) REFERENCES status(id)
 );
