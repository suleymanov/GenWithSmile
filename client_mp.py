import os
import sys
from json import loads
from multiprocessing import Process, Pool, TimeoutError

from comb_lib_functions import get_combinations


def _process2(params):
	i = params[0]
	mol_smiles_list = params[1]
	src_atom = params[2]
	trg_atom = params[3]
	result_fn = params[4]
	for j in xrange(i, len(mol_smiles_list)):
		# print('i = {},\tj = {}'.format(i, j))
		_combinations = get_combinations(mol_smiles_list[i], mol_smiles_list[j], src_atom, trg_atom)
		with open(result_fn.format(i, j), 'w') as f:
			f.write('\n'.join(_combinations))


def main(args):
	settings_fn = args[0]
	beg_ind = int(args[1])
	end_ind = int(args[2])
	num_proc = int(args[3])
	with open(settings_fn) as f:
		settings = loads(f.read())
	smiles_fn = settings['smiles_file_name']
	source_atom = settings['source_atom']
	target_atom = settings['target_atom']
	result_fn = settings['result_file_name']
	with open(smiles_fn) as f:
		mol_smiles_list = f.read().splitlines()
	
	def _process(i):
		for j in xrange(i, len(mol_smiles_list)):
			_combinations = get_combinations(mol_smiles_list[i], mol_smiles_list[j], src_atom, trg_atom)
			with open(result_fn.format(i, j), 'w') as f:
				f.write('\n'.join(_combinations))

	pack_i = lambda i: (i, mol_smiles_list, source_atom, target_atom, result_fn)

	for i_start in xrange(beg_ind, end_ind, num_proc):
		pool = Pool(processes=num_proc)
		# pool.map(_process, range(i_start, i_start+num_proc))
		pool.map(_process2, map(pack_i, xrange(i_start, i_start+num_proc)))


if __name__ == '__main__':
	main(sys.argv[1:])
