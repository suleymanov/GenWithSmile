import os
import unittest
import tempfile
from random import sample

from gws.src import SourceTargetClient
from gws.src.core import get_unique_mols, MoleculeHandler
from gws.src.io import IOUtils

from utils import ClientBaseTestCase
from utils import clean_files, clean_paths


class ClientsSourceTargetTests(ClientBaseTestCase):
	def test_1(self):
		print('Test 1: atoms as inputs')
		self.assertTrue(self._run_test_case(self.atoms, self.atoms))
		self.assertTrue(self._run_test_case(self.atoms, self.basic, attach=False, merge=True))
		self.assertTrue(self._run_test_case(self.atoms, self.groups))
		self.assertTrue(self._run_test_case(self.atoms, self.rings, attach=False, merge=True))

	def test_2(self):
		print('Test 2: basic structs as inputs')
		self.assertTrue(self._run_test_case(self.basic, self.atoms))
		self.assertTrue(self._run_test_case(self.basic, self.basic, attach=False, merge=True))
		self.assertTrue(self._run_test_case(self.basic, self.groups))
		self.assertTrue(self._run_test_case(self.basic, self.rings, attach=False, merge=True))

	def test_3(self):
		print('Test 3: groups as inputs')
		self.assertTrue(self._run_test_case(self.groups, self.atoms))
		self.assertTrue(self._run_test_case(self.groups, self.basic, attach=False, merge=True))
		self.assertTrue(self._run_test_case(self.groups, self.groups))
		self.assertTrue(self._run_test_case(self.groups, self.rings, attach=False, merge=True))

	def test_4(self):
		print('Test 4: rings as inputs')
		self.assertTrue(self._run_test_case(self.rings, self.atoms))
		self.assertTrue(self._run_test_case(self.rings, self.basic, attach=False, merge=True))
		self.assertTrue(self._run_test_case(self.rings, self.groups))
		self.assertTrue(self._run_test_case(self.rings, self.rings, attach=False, merge=True))

	def _run_test_case(self, structures, addons, num_iter=1, num_add=1, two_point=False,
		attach=True, merge=False):
		status = True
		num_tests = min(self.num_repeats, len(structures))
		sources = sample(structures, num_tests)
		targets = sample(structures, num_tests)
		for src, trg in zip(sources, targets):
			source = {'smiles': src}
			target = {'smiles': trg}
			iterations = []
			for _ in xrange(num_iter):
				adds = map(lambda s: {'smiles': s}, sample(addons, num_add))
				iterations.append({
					'attach': {'one_point': attach, 'two_point': two_point, 'addons': adds},
					'merge': {'one_point': merge, 'two_point': two_point, 'addons': adds}})
			config_fn = create_config(source, target, iterations)
			client = SourceTargetClient.from_file(config_fn)
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
			# handlers = [MoleculeHandler(s) for s in mols]
			status = status and len(handlers) == len(get_unique_mols(map(lambda x: x.mol.graph, handlers)))

		return status


def create_config(source, target, iterations):
	temp_fn = tempfile.mkstemp()[1]
	temp_pn = tempfile.mkdtemp()
	json_config = {
		'source': source, 'target': target, 'iterations': iterations, 'output': {'path': temp_pn}}
	IOUtils.write_json(temp_fn, json_config)
	return temp_fn


if __name__ == '__main__':
	unittest.main()
