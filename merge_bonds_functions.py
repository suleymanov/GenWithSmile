from rdkit import Chem

from utils import (
	rdkitmol2graph, get_nonisomorphic_positions_2, is_isomorph, is_isomorph_gk, vectorize_mol_graphs
)
from comb_lib_functions import get_unique_mols
from combine import merge_bond


def get_unique_inds(mol_smiles):
	"""
	Get positions of bonds which provide unique result under some sort of modification.
	:param mol_smiles: str
	:return: list of list of int
	"""
	mol = Chem.MolFromSmiles(mol_smiles)
	if not mol:
		raise StandardError('Cannot convert SMILES to rdkit molecule.')
	if mol.GetNumAtoms() == 1:
		return None
	mol_graph = rdkitmol2graph(mol)
	bonds_positions = map(
		lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), 
		filter(lambda b: b.IsInRing(), mol.GetBonds())
	)
	return [bonds_positions[ind] 
		for ind in get_nonisomorphic_positions_2(mol_graph, bonds_positions)]


def get_combinations(mol_smiles_1, mol_smiles_2):
	"""
	Generate unique molecules in SMILES format by pairwise bonds merging.
	:param mol_smiles_1: str
	:param mol_smiles_2: str
	:return: list of str
	"""
	unique_inds_1 = get_unique_inds(mol_smiles_1)
	unique_inds_2 = get_unique_inds(mol_smiles_2)
	results = []
	for bond_pos_1 in unique_inds_1:
		for bond_pos_2 in unique_inds_2:
			combinations = merge_bond(mol_smiles_1, mol_smiles_2, bond_pos_1, bond_pos_2)
			results.extend(combinations)
			results = get_unique_mols(results)
	return results
