import os
import sys
from json import loads

from comb_lib_functions import get_unique_inds, get_unique_mols, get_combinations


def main(args):
	settings_fn = args[0]
	with open(settings_fn) as f:
		settings = loads(f.read())
	smiles_fn = settings['smiles_file_name']
	source_atom = settings['source_atom']
	target_atom = settings['target_atom']
	result_fn = settings['result_file_name']
	with open(smiles_fn) as f:
		mol_smiles_list = f.read().splitlines()
	# for i in xrange(len(mol_smiles_list)):
        # for i in xrange(500, len(mol_smiles_list)):
	for i in xrange(125, 150):
		for j in xrange(i, len(mol_smiles_list)):
			print('i = {},\tj = {}'.format(i, j))
			_combinations = get_combinations(
				mol_smiles_list[i], mol_smiles_list[j], source_atom, target_atom)
			with open(result_fn.format(i, j), 'w') as f:
				f.write('\n'.join(_combinations))


if __name__ == '__main__':
	main(sys.argv[1:])
