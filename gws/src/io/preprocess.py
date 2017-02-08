from __future__ import absolute_import

import os
import sys
import tempfile
from math import sqrt

from rdkit import Chem

from ..core import get_unique_mols
from ..core.algs import get_unique_coords, rdkitmol2graph
from .utils import IOUtils


def basic(args):
	"""
	Basic input files preprocessing:
	- convert *.smi to *.json along with creating appropriate fields;
	- select 'unique' attach/merge positions (possibly filtering through provided atoms symbols);
	- select all unique molecules.
	:param args: list of str
	:return: None
	"""
	smiles_fn_in = args[0]
	json_fn_out = args[1]
	atoms = args[2].split(',') if len(args) == 3 else []

	smiles_list = IOUtils.read_smi(smiles_fn_in)
	mol_graph_list = map(
		rdkitmol2graph, 
		filter(lambda x: x is not None, map(Chem.MolFromSmiles, smiles_list)))
	unique_inds = get_unique_mols(mol_graph_list)

	processed = []
	for ind in unique_inds:
		mol_graph = mol_graph_list[ind]
		pos = sorted(map(lambda x: x[0], get_unique_coords(
			mol_graph, 1, map(lambda x: [x], xrange(mol_graph.number_of_nodes())))))
		if atoms:
			pos = filter(lambda node_key: mol_graph.node[node_key]['label'] in atoms, mol_graph.nodes())
		mol_data = {'smiles': smiles_list[ind], 'attach_pos': pos, 'merge_pos': pos}
		processed.append(mol_data)

	IOUtils.write_json(json_fn_out, processed)


def select_3d(args):
	"""
	3D representation preprocessing:
	- convert *.mol2 to *.json with creating appropriate fields;
	- select atom indexes that fall into specified spheres;
	- select 'unique' attach/merge positions (possibly filtering through provided atoms symbols);
	- select all unique molecules.
	:param args: list of str
	:return: None
	"""
	def split_multiple_mol2(fn):
		data = IOUtils.read_file(fn).splitlines()
		start_inds = map(lambda x: x[0], filter(lambda x: '<TRIPOS>MOLECULE' in x[1], enumerate(data)))
		end_inds = start_inds[1:][::] + [len(data)]
		return map(lambda coord: data[coord[0]:coord[1]], zip(start_inds, end_inds))

	def process_single_mol2(content, spheres_settings, atoms):
		def get_line(lines, patt):
			line_ind = filter(lambda i: patt in lines[i], xrange(len(lines)))
			assert len(line_ind) == 1
			return line_ind[0]

		fn = tempfile.mkstemp()[1]
		IOUtils.write_file(fn, '\n'.join(content) + '\n')
		calc_dist = lambda coord, coord2: sqrt(sum((x - y) * (x - y) for x, y in zip(coord, coord2)))
		smiles = Chem.MolToSmiles(Chem.MolFromMol2File(fn))
		mol_graph = rdkitmol2graph(Chem.MolFromSmiles(smiles))
		os.remove(fn)
		
		start, end = get_line(content, '<TRIPOS>ATOM'), get_line(content, '<TRIPOS>BOND')
		atoms_data = map(lambda x: x.split(), content[(start+1):end])
		atoms_data = filter(lambda x: x[5].split('.')[0] != 'H', atoms_data)
		pos = []
		for i, item in enumerate(atoms_data):
			if atoms and item[5] not in atoms:
				continue
			coords = map(float, item[2:5])
			if any(map(lambda x: calc_dist(x['center'], coords) <= x['dist'], spheres_settings)):
				pos.append(i)
		
		pos = sorted(map(lambda x: x[0], get_unique_coords(mol_graph, 1, map(lambda x: [x], pos))))
		return {'smiles': smiles, 'attach_pos': pos, 'merge_pos': pos}, mol_graph

	mol2_fn_in = args[0]
	json_fn_out = args[1]
	spheres_settings = IOUtils.read_json(args[2])
	atoms = args[3].split(',') if len(args) == 4 else []

	assert len(spheres_settings) > 0
	assert all(map(lambda x: x['dist'] >= 0.0, spheres_settings))
	assert all(map(lambda x: len(x['center']) == 3, spheres_settings))

	splitted_data = split_multiple_mol2(mol2_fn_in)
	processed = []
	for content in splitted_data:
		processed.append(process_single_mol2(content, spheres_settings, atoms))
	IOUtils.write_json(
		json_fn_out,
		map(lambda ind: processed[ind][0], get_unique_mols(map(lambda x: x[1], processed))))


def main(args):
	option = args[0]
	if option == '--basic':
		basic(args[1:])
	if option == '--select3d':
		select_3d(args[1:])


if __name__ == '__main__':
	main(sys.argv[1:])
