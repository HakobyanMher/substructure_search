import pytest
from src.search import substructure_search


def test_benzene_example():
    data = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    sub = "c1ccccc1"
    out = substructure_search(data, sub)
    assert out == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]

def test_invalid_substructure_raises():
    data = ["CCO", "c1ccccc1"]
    with pytest.raises(ValueError):
        substructure_search(data, "not_a_smiles")

def test_invalid_molecule_is_skipped():
    data = ["CCO", "###not-smiles###", "c1ccccc1"]
    sub = "c1ccccc1"
    out = substructure_search(data, sub)
    assert out == ["c1ccccc1"]

def test_empty_input():
    assert substructure_search([], "C") == []
