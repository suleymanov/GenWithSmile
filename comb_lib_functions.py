from itertools import izip

from rdkit import Chem

from utils import rdkitmol2graph
from utils import get_nonisomorphic_positions
from utils import is_isomorph
from combine import combine_mols


def get_unique_inds(mol_smiles):
	"""
	Get positions of atoms which provide unique result under some sort of modification.
	:param mol_smiles: str
	:return: list of int
	"""
	mol = Chem.MolFromSmiles(mol_smiles)
	if mol.GetNumAtoms() == 1:
		return [0]
	mol_graph = rdkitmol2graph(mol)
	positions = list(xrange(len(mol.GetAtoms())))
	return get_nonisomorphic_positions(mol_graph, positions)


def get_unique_mols(mol_smiles_list, use_gk=False):
	"""
	Filter non-unique molecules.
	:param mol_smiles_list: list of str
	:return: list of str
	"""
	if not mol_smiles_list:
		return []
	unique_inds = [0]
	mol_graphs = map(rdkitmol2graph, map(Chem.MolFromSmiles, mol_smiles_list))
	for i, mol_graph in enumerate(mol_graphs[1:]):
		if any(is_isomorph(mol_graph, mol_graphs[ind]) for ind in unique_inds):
			continue
		unique_inds.append(i + 1)
	return [mol_smiles_list[ind] for ind in unique_inds]


def get_combinations2(args):
	return args


def get_combinations(mol_smiles_1, mol_smiles_2, src_atom, trg_atom):
	"""
	Generate unique and chemically valid combiinations of two molecules.
	:param mol_smiles_1: str
	:param mol_smiles_2: str
	:param src_atom: str
	:param trg_atom: str
	:return: list of str
	"""
	unique_inds_1 = get_unique_inds(mol_smiles_1)
	unique_inds_2 = get_unique_inds(mol_smiles_2)
	results = []
	for ind_1 in unique_inds_1:
		for ind_2 in unique_inds_2:
			combinations = combine_mols(
				mol_smiles_1, mol_smiles_2, ind_1, ind_2, src_atom, trg_atom)
			results.extend(combinations)
			results = get_unique_mols(results)
	return results
