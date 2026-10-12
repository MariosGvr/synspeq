"""Molecular identity policy.

Identity = RDKit canonical isomeric SMILES of the molecule as read (stereochemistry kept, hydrogens
removed, no desalting, no neutralisation, no tautomer normalisation, no largest-fragment selection).
Nothing is cleaned up silently: a record that RDKit cannot sanitise is reported as invalid, not repaired.
Conformers of one structure therefore share one identity.
"""
from __future__ import annotations
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
POLICY_ID = "rdkit-canonical-isomeric-v1"


def identity_from_mol(mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def identity_from_smiles(smiles: str):
    """Return (identity, reason). identity is None when the string cannot be parsed/sanitised."""
    if smiles is None or not str(smiles).strip():
        return None, "empty"
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None, "unparsable_smiles"
    if mol.GetNumAtoms() == 0:
        return None, "no_atoms"
    return identity_from_mol(mol), None


def n_fragments(identity: str) -> int:
    return identity.count(".") + 1
