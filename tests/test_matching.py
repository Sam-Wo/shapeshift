"""Structure-matching unit tests (CPU, seconds). Covers the cases that drove the design of overlap.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from rdkit import Chem  # noqa: E402

from overlap import METALS, graph_hash, norm_name, std_key, std_mol  # noqa: E402


def formula(smiles):
    return Chem.MolToInchi(std_mol(smiles)).split("/")[1]


def test_salt_and_charge_are_stripped():
    free = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"   # imatinib
    salt = free + ".CS(=O)(=O)O"                                                      # imatinib mesylate
    assert std_key(free) == std_key(salt)


def test_stereo_is_flattened():
    # JUMP InChIKeys are stereo-flat, so enantiomers must collapse to one key
    assert std_key("C[C@H](N)C(=O)O") == std_key("C[C@@H](N)C(=O)O") == std_key("CC(N)C(=O)O")


def test_shifted_double_bond_matches_on_graph_not_key():
    # dexamethasone as Tahoe draws it vs the delta-5 diene JUMP stores for positive control JCP2022_025848
    tahoe = "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C"
    jump = "CC1CC2C3CC=C4CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO"
    assert std_key(tahoe) != std_key(jump)              # InChIKey alone misses it ...
    assert formula(tahoe) == formula(jump)              # ... same formula
    assert graph_hash(std_mol(tahoe)) == graph_hash(std_mol(jump))   # ... same heavy-atom graph -> tier 2b


def test_regioisomer_analog_is_not_matched():
    # belzutifan (2,3-difluoro) vs the 2,2-difluoro analog in JUMP: same formula, different graph
    belzutifan = "CS(=O)(=O)C1=C2C(C(C(C2=C(C=C1)OC3=CC(=CC(=C3)C#N)F)F)F)O"
    analog = "CS(=O)(=O)c1ccc(Oc2cc(F)cc(C#N)c2)c2c1C(O)C(F)(F)C2"
    assert formula(belzutifan) == formula(analog)
    assert graph_hash(std_mol(belzutifan)) != graph_hash(std_mol(analog))


def test_metal_complex_is_flagged():
    auranofin = "CCP(CC)CC.CC(=O)OCC1C(C(C(C(O1)[S-])OC(=O)C)OC(=O)C)OC(=O)C.[Au+]"
    assert any(a.GetSymbol() in METALS for a in Chem.MolFromSmiles(auranofin).GetAtoms())


def test_name_normalisation():
    assert norm_name("Erdafitinib ") == norm_name("erdafitinib")
    assert norm_name("Fluvoxamine (maleate)") == norm_name("fluvoxamine")
