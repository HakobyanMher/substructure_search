from rdkit import Chem


def substructure_search(molecules: list[str], substructure: str) -> list[str]:
    query = Chem.MolFromSmiles(substructure)
    if query is None:
        raise ValueError(f"Invalid substructure SMILES: {substructure!r}")

    hits: list[str] = []
    for smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        if mol.HasSubstructMatch(query):
            hits.append(smiles)
    return hits
