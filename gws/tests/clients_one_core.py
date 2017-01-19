import os
import unittest
import tempfile
from random import sample

from gws.src import OneCoreClient
from gws.src.core import get_unique_mols, MoleculeHandler
from gws.src.io import IOUtils

from utils import ClientBaseTestCase
from utils import clean_files, clean_paths


class ClientsOneCoreTests(ClientBaseTestCase):
	def test_1(self):
		print('Test 1: atoms as cores')
		self.assertTrue(self._run_test_case(self.atoms, self.atoms))
		self.assertTrue(self._run_test_case(self.atoms, self.basic))
		self.assertTrue(self._run_test_case(self.atoms, self.groups))
		self.assertTrue(self._run_test_case(self.atoms, self.rings))

	def test_2(self):
		print('Test 2: basic structs as cores')
		self.assertTrue(self._run_test_case(self.basic, self.atoms))
		self.assertTrue(self._run_test_case(self.basic, self.basic))
		self.assertTrue(self._run_test_case(self.basic, self.groups))
		self.assertTrue(self._run_test_case(self.basic, self.rings))

	def test_3(self):
		print('Test 3: groups as cores')
		self.assertTrue(self._run_test_case(self.groups, self.atoms))
		self.assertTrue(self._run_test_case(self.groups, self.basic))
		self.assertTrue(self._run_test_case(self.groups, self.groups))
		self.assertTrue(self._run_test_case(self.groups, self.rings))

	def test_4(self):
		print('Test 4: rings as cores')
		self.assertTrue(self._run_test_case(self.rings, self.atoms))
		self.assertTrue(self._run_test_case(self.rings, self.basic))
		self.assertTrue(self._run_test_case(self.rings, self.groups))
		self.assertTrue(self._run_test_case(self.rings, self.rings))

	def test_5(self):
		print('Test 5: two iterations')
		self.assertTrue(self._run_test_case(self.basic, self.groups, num_iter=2))
		self.assertTrue(self._run_test_case(self.groups, self.rings, num_iter=2))
		self.assertTrue(self._run_test_case(self.rings, self.basic, num_iter=2))

	def test_6(self):
		print('Test 6: two modifications')
		self.assertTrue(self._run_test_case(self.basic, self.groups, max=2))
		self.assertTrue(self._run_test_case(self.groups, self.rings, max=2))
		self.assertTrue(self._run_test_case(self.rings, self.basic, max=2))

	def test_7(self):
		print('Test 7: two addons')
		self.assertTrue(self._run_test_case(self.basic, self.groups, num_add=2))
		self.assertTrue(self._run_test_case(self.groups, self.rings, num_add=2))
		self.assertTrue(self._run_test_case(self.rings, self.basic, num_add=2))

	def test_8(self):
		print('Test 7: two-point')
		self.assertTrue(self._run_test_case(self.basic, self.groups, two_point=True))
		self.assertTrue(self._run_test_case(self.groups, self.rings, two_point=True))
		self.assertTrue(self._run_test_case(self.rings, self.basic, two_point=True))

	def _run_test_case(self, structures, addons, num_iter=1, max=1, num_add=1, two_point=False):
		status = True
		num_tests = min(self.num_repeats, len(structures))
		cores = sample(structures, num_tests)
		for c in cores:
			core = {'smiles': c}
			iterations = []
			for _ in xrange(num_iter):
				adds = map(lambda s: {'smiles': s}, sample(addons, num_add))
				iterations.append({
					'attach': {'two_point': two_point, 'addons': adds},
					'merge': {'two_point': two_point, 'addons': adds},
					'max': max})
			config_fn = create_config(core, iterations)
			client = OneCoreClient.from_file(config_fn)
			self.files_to_remove.append(config_fn)
			self.paths_to_remove.append(client._config.output.path)
			client.process()
			status = status and self._all_unique(client._config)
		return status

	def _all_unique(self, run_config):
		status = True
		path = run_config.output.path

		for i in xrange(len(run_config.iterations)):
			mols = reduce(
				lambda res, item: res + IOUtils.read_smi(path + os.sep + item),
				filter(lambda x: 'iter_{}'.format(i + 1) in x, os.listdir(path)), [])
			handlers = map(MoleculeHandler, mols)
			status = status and len(handlers) == len(get_unique_mols(map(lambda x: x.mol.graph, handlers)))

		return status


def create_config(core, iterations):
	temp_fn = tempfile.mkstemp()[1]
	temp_pn = tempfile.mkdtemp()
	json_config = {'core': core, 'iterations': iterations, 'output': {'path': temp_pn}}
	IOUtils.write_json(temp_fn, json_config)
	return temp_fn


if __name__ == '__main__':
	unittest.main()
