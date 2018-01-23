CREATE TABLE status (
    id TEXT NOT NULL UNIQUE,
    status INTEGER NOT NULL
);

CREATE TABLE result (
    id TEXT NOT NULL UNIQUE,
    chem_index INTEGER NOT NULL,
    smiles TEXT NOT NULL,
    hydrolyzability INTEGER
 );
