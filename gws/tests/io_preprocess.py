import os
import unittest
import tempfile
from math import sqrt

from rdkit import Chem

from gws.src.io import IOUtils
from gws.src.io import preprocess
from gws.src.core.algs import get_unique_coords, rdkitmol2graph


class PreprocessSelect3DTests(unittest.TestCase):
	def setUp(self):
		self.settings_fn = tempfile.mkstemp()[1]
		self.results_fn = tempfile.mkstemp()[1]
		self.mol2_fn = 'gws/tests/resources/3d_example.mol2'
		mol_graph = rdkitmol2graph(Chem.MolFromMol2File(self.mol2_fn))
		self.unique_coords = map(lambda x: x[0], get_unique_coords(
			mol_graph, 1, map(lambda x: [x], xrange(mol_graph.number_of_nodes()))))

		contents = IOUtils.read_file(self.mol2_fn).splitlines()
		i_start = map(
			lambda x: x[0], 
			filter(lambda x: '<TRIPOS>ATOM' in x[1], enumerate(contents)))[0] + 1
		i_end = map(
			lambda x: x[0], 
			filter(lambda x: '<TRIPOS>BOND' in x[1], enumerate(contents)))[0]
		self.atom_data = filter(
			lambda strs: strs[5].split('.')[0] != 'H', 
			map(lambda x: x.split(), contents[i_start:i_end]))
		self.max_dist = max(sqrt(
			sum((x - y) * (x - y) 
				for x, y in zip(map(float, atom_item[2:5]), map(float, atom_item2[2:5]))))
			for atom_item in self.atom_data for atom_item2 in self.atom_data)

	def tearDown(self):
		os.remove(self.settings_fn)
		os.remove(self.results_fn)

	def test1(self):
		for i, atom_item in enumerate(self.atom_data):
			settings = [{'center': map(float, atom_item[2:5]), 'dist': 0.0}]
			IOUtils.write_json(self.settings_fn, settings)
			preprocess.select_3d([self.mol2_fn, self.results_fn, self.settings_fn])
			result = IOUtils.read_json(self.results_fn)[0]
			self.assertTrue(len(result['attach_pos']) == 1)
			self.assertTrue(result['attach_pos'][0] == i)

	def test2(self):
		settings = [{'center': map(float, self.atom_data[0][2:5]), 'dist': self.max_dist}]
		IOUtils.write_json(self.settings_fn, settings)
		preprocess.select_3d([self.mol2_fn, self.results_fn, self.settings_fn])
		self.assertTrue(
			len(IOUtils.read_json(self.results_fn)[0]['attach_pos']) == len(self.unique_coords))

	def test3(self):
		res = []
		for i in xrange(10):
			settings = [{
				'center': map(float, self.atom_data[0][2:5]), 
				'dist': (i + 1) * self.max_dist / 10}]
			IOUtils.write_json(self.settings_fn, settings)
			preprocess.select_3d([self.mol2_fn, self.results_fn, self.settings_fn])
			res.append(len(IOUtils.read_json(self.results_fn)[0]['attach_pos']))
		self.assertTrue(all(x <= y for x, y in zip(res[:-1], res[1:])))


if __name__ == '__main__':
	unittest.main()
