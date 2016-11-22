import os
import sys
from json import loads

from comb_lib_functions import get_combinations

# from rdkit import Chem


def main(args):
	settings_fn = args[0]
	beg_ind = int(args[1])
	end_ind = int(args[2])
	with open(settings_fn) as f:
		settings = loads(f.read())
	smiles_fn = settings['smiles_file_name']
	source_atom = settings['source_atom']
	target_atom = settings['target_atom']
	result_fn = settings['result_file_name']
	with open(smiles_fn) as f:
		mol_smiles_list = f.read().splitlines()
	for i in xrange(beg_ind, end_ind):
		for j in xrange(i, len(mol_smiles_list)):
			print('i = {},\tj = {}'.format(i, j))
			_combinations = get_combinations(
				mol_smiles_list[i], mol_smiles_list[j], source_atom, target_atom)
			with open(result_fn.format(i, j), 'w') as f:
				f.write('\n'.join(_combinations) + '\n')


if __name__ == '__main__':
	main(sys.argv[1:])
