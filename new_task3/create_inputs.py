import os
import sys
from itertools import izip
from json import dumps, loads

from rdkit import Chem


def _create_inputs_1(mols_data):
	"""
	Single rings.
	"""
	inputs = []
	for item in filter(lambda x: x['single_ring'], mols_data):
		mol = Chem.MolFromSmiles(item['smiles'])
		assert all(map(lambda x: mol.GetAtomWithIdx(x).IsInRing(), item['pos']))
		inputs.append({'smiles': item['smiles'], 'attach_pos': sorted(item['pos']), 'atoms': ['C']})
	return inputs


def _create_inputs_2(mols_data):
	"""
	All rings.
	"""
	inputs = []
	for item in mols_data:
		mol = Chem.MolFromSmiles(item['smiles'])
		assert all(map(lambda x: mol.GetAtomWithIdx(x).IsInRing(), item['pos']))
		inputs.append({'smiles': item['smiles'], 'merge_pos': sorted(item['pos']), 'atoms': ['C']})
	return inputs


def main(args):
	mols_fn = args[0]  # molecules in SMILES format
	fn_1 = args[1]
	fn_2 = args[2]
	with open(mols_fn) as f:
		mols_data = map(
			lambda strs: {
				'smiles': strs[0], 
				'pos': sorted(map(int, strs[2].split(','))), 
				'single_ring': len(strs) == 4}, 
			map(lambda line: line.split(), f.read().splitlines()))
	list_1 = _create_inputs_1(mols_data)
	list_2 = _create_inputs_2(mols_data)
	with open(fn_1, 'w') as f:
		f.write(dumps(list_1, indent=4))
	with open(fn_2, 'w') as f:
		f.write(dumps(list_2, indent=4))


if __name__ == '__main__':
	main(sys.argv[1:])
