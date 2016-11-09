import os
import sys

from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from functions import get_combinations, get_unique_mols


def main(args):
	"""
	Main demo script.
	:param args: list of str
	:return: None
	"""
	mol_smiles_1 = args[0]
	mol_smiles_2 = args[1]
	num_iter = int(args[2])
	out_fn_template = args[3]
	out_fn_results = args[4]
	mol_smiles_list = _generate(mol_smiles_1, mol_smiles_2, num_iter)
	_draw_mols(mol_smiles_list, out_fn_template)
	with open(out_fn_results, 'w') as f:
		f.write('\n'.join(mol_smiles_list) + '\n')


def _generate(mol_smiles_1, mol_smiles_2, num_iter):
	hosts = [mol_smiles_1]
	processed_hosts = []
	addons = [mol_smiles_2]
	results = []
	for _ in xrange(num_iter):
		while len(hosts) > 0:
			host_smiles = hosts.pop()
			processed_hosts.append(host_smiles)
			for addon_smiles in addons:
				_combinations = get_combinations(host_smiles, addon_smiles)
				results.extend(_combinations)
				results = get_unique_mols(results)
		hosts = filter(lambda x: x not in processed_hosts, results)
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
