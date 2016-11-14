import os
import sys

from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from comb_lib_functions import get_unique_mols
from merge_bonds_functions import get_combinations


def main(args):
	"""
	Main demo script.
	:param args: list of str
	:return: None
	"""
	mol_smiles_fn = args[0]
	out_fn_template = args[1]
	out_fn_results = args[2]
	num_iter = int(args[3])
	with open(mol_smiles_fn) as f:
		mol_smiles_list = f.read().splitlines()
	results = _generate(mol_smiles_list, num_iter)
	_draw_mols(results, out_fn_template)
	with open(out_fn_results, 'w') as f:
		f.write('\n'.join(results) + '\n')


def _generate_1_iter(mol_smiles_list):
	results = []
	for mol_smiles_1 in mol_smiles_list:
		for mol_smiles_2 in mol_smiles_list:
			_combinations = get_combinations(mol_smiles_1, mol_smiles_2)
			results.extend(_combinations)
			results = get_unique_mols(results)
	return results


def _generate(mol_smiles_list, num_iter):
	results = []
	hosts = mol_smiles_list[::]
	for _ in xrange(num_iter):
		for mol_smiles_1 in hosts:
			for mol_smiles_2 in hosts:
				_combinations = get_combinations(mol_smiles_1, mol_smiles_2)
				results.extend(_combinations)
				results = get_unique_mols(results, use_gk=True)
		hosts.extend(results)
		hosts = list(set(hosts))
	return results


def _draw_mols(mol_smiles_list, out_fn_template):
	num_rows = 3
	mols_per_row = 3

	rdkit_mols = map(Chem.MolFromSmiles, mol_smiles_list)

	for i, i_start in enumerate(xrange(0, len(rdkit_mols), num_rows * mols_per_row)):
		subset = rdkit_mols[i_start:i_start + num_rows * mols_per_row]
		fn = '{}_{}.png'.format(out_fn_template, i + 1)
		img = Draw.MolsToGridImage(subset, molsPerRow=mols_per_row, subImgSize=(100, 100))
		img.save(fn)


if __name__ == '__main__':
	main(sys.argv[1:])
