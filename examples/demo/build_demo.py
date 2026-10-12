"""Builds the small demo input in this folder: molecules.sdf, routes/*.json and building_blocks.smi.
The planner records are written by hand in AiZynthFinder's JSON layout; they illustrate the input format and are not
planner output. Run:  python examples/demo/build_demo.py"""
import json, os
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "routes"), exist_ok=True)

MOLS = [  # (name, smiles, pocket)
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "pocket_A"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", "pocket_A"),
    ("benzocaine", "CCOC(=O)c1ccc(N)cc1", "pocket_A"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "pocket_B"),
    ("ethanol", "CCO", "pocket_B"),
    ("ethanol_again", "OCC", "pocket_B"),            # the same molecule: counted as a duplicate
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "pocket_B"),   # no route in the demo records
]


def mol(smi, stock=False, children=None):
    n = {"type": "mol", "smiles": smi, "in_stock": stock}
    if children is not None:
        n["children"] = children
    return n


def rxn(*kids):
    return {"type": "reaction", "children": list(kids)}


def record(name, smi, tree=None, steps=0):
    return {"mol_id": name, "smiles": smi, "solvable": tree is not None, "n_steps": steps if tree else None,
            "trees": [tree] if tree else [], "error": None}


R = {
    "aspirin": record("aspirin", "CC(=O)Oc1ccccc1C(=O)O",
                      mol("CC(=O)Oc1ccccc1C(=O)O", children=[rxn(mol("O=C(O)c1ccccc1O", True), mol("CC(=O)OC(C)=O", True))]), 1),
    "paracetamol": record("paracetamol", "CC(=O)Nc1ccc(O)cc1",
                          mol("CC(=O)Nc1ccc(O)cc1", children=[rxn(mol("Nc1ccc(O)cc1", True), mol("CC(=O)OC(C)=O", True))]), 1),
    "benzocaine": record("benzocaine", "CCOC(=O)c1ccc(N)cc1",
                         mol("CCOC(=O)c1ccc(N)cc1", children=[rxn(mol("Nc1ccc(C(=O)O)cc1", True), mol("CCO", True))]), 1),
    "ibuprofen": record("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
                        mol("CC(C)Cc1ccc(C(C)C(=O)O)cc1", children=[rxn(
                            mol("CC(C)Cc1ccc(C(C)C#N)cc1", children=[rxn(mol("CC(C)Cc1ccc(C(C)Cl)cc1", True), mol("N#C[Na]", True))]))]), 2),
    "ethanol": record("ethanol", "CCO", mol("CCO", True), 0),
    "ethanol_again": record("ethanol_again", "CCO", mol("CCO", True), 0),
    "caffeine": record("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
}
for k, v in R.items():
    json.dump(v, open(os.path.join(HERE, "routes", f"{k}.json"), "w"), indent=1)

blocks = ""
for name, smi, pocket in MOLS:
    m = Chem.MolFromSmiles(smi)
    m.SetProp("_Name", name)
    blocks += Chem.MolToMolBlock(m) + f">  <pocket_id>\n{pocket}\n\n$$$$\n"
bad = Chem.MolFromSmiles("C[N](C)(C)(C)C", sanitize=False)          # pentavalent nitrogen: reported as invalid
bad.SetProp("_Name", "bad_record")
blocks += Chem.MolToMolBlock(bad, kekulize=False) + ">  <pocket_id>\npocket_B\n\n$$$$\n"
open(os.path.join(HERE, "molecules.sdf"), "w").write(blocks)

# a tiny purchasing catalogue: the nitrile intermediate of the ibuprofen route is stocked
open(os.path.join(HERE, "building_blocks.smi"), "w").write("CC(C)Cc1ccc(C(C)C#N)cc1\nNc1ccc(C(=O)O)cc1\n")
print("demo written to", os.path.relpath(HERE))
