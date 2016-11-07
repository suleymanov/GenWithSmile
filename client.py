import os
import sys

from functions import get_unique_inds, get_unique_mols, get_combinations


def _unique_inds(args):
	# Some examples:

	# C1CCCC1(N)
	# C1=C(C=CC=C1)C2=CC=CC=C2
	# C1CCCCC1C3CCCC(C2CCCCC2)C3
	# C1CCCCC1C3CCCC(=C2CCCCC2)C3
	mol_smiles = args[0]
	return get_unique_inds(mol_smiles)


def _unique_mols(args):
	fn_in = args[0]
	fn_out = args[1]
	use_gk = args[2] == '--use_gk'
	with open(fn_in) as f:
		mol_smiles_list = f.read().splitlines()

	mol_smiles_list_unique = get_unique_mols(mol_smiles_list, use_gk)
	with open(fn_out, 'w') as f:
		f.write('\n'.join(mol_smiles_list_unique) + '\n')


def _combinations(args):
	mol_smiles_1 = args[0]
	mol_smiles_2 = args[1]
	fn_out = args[2]
	combinations = get_combinations(mol_smiles_1, mol_smiles_2)
	with open(fn_out, 'w') as f:
		f.write('\n'.join(combinations) + '\n')


def main(args):
	option = args[0]
	if option == '--unique_inds':
		print(_unique_inds(args[1:]))
		return
	if option == '--unique_mols':
		print(_unique_mols(args[1:]))
		return
	if option == '--combinations':
		print(_combinations(args[1:]))
		return
	raise KeyError('Unknown option.')


if __name__ == '__main__':
	main(sys.argv[1:])
