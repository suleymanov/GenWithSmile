import os
import unittest
import tempfile
from random import sample

from gws.src import CombinationsClient
from gws.src.core import get_unique_mols, MoleculeHandler
from gws.src.io import IOUtils

from utils import ClientBaseTestCase
from utils import clean_files, clean_paths


class ClientsCombinationsTests(ClientBaseTestCase):
	def test_1(self):
		print('Test 1: X * X')
		self.assertTrue(self._run_test_case(self.atoms))
		self.assertTrue(self._run_test_case(self.basic))
		self.assertTrue(self._run_test_case(self.groups))
		self.assertTrue(self._run_test_case(self.rings))

	def test_2(self):
		print('Test 2: X * Y')
		self.assertTrue(self._run_test_case(self.atoms, self.basic))
		self.assertTrue(self._run_test_case(self.basic, self.groups))
		self.assertTrue(self._run_test_case(self.groups, self.rings))
		self.assertTrue(self._run_test_case(self.rings, self.atoms))

	def test_3(self):
		print('Test 3: X * L * X')
		self.assertTrue(self._run_test_case(self.atoms, linkers=self.basic))
		self.assertTrue(self._run_test_case(self.basic, linkers=self.groups))
		self.assertTrue(self._run_test_case(self.groups, linkers=self.rings))
		self.assertTrue(self._run_test_case(self.rings, linkers=self.atoms))

	def test_4(self):
		print('Test 4: X * L * Y')
		self.assertTrue(self._run_test_case(self.atoms, self.basic, self.groups))
		self.assertTrue(self._run_test_case(self.basic, self.groups, self.rings))
		self.assertTrue(self._run_test_case(self.groups, self.rings, self.atoms))
		self.assertTrue(self._run_test_case(self.rings, self.atoms, self.basic))

	def test_5(self):
		print('Test 5: X * L * Y (two-point)')
		self.assertTrue(self._run_test_case(self.atoms, self.basic, self.groups, one_point=False, two_point=True))
		self.assertTrue(self._run_test_case(self.basic, self.groups, self.rings, one_point=False, two_point=True))
		self.assertTrue(self._run_test_case(self.groups, self.rings, self.atoms, one_point=False, two_point=True))
		self.assertTrue(self._run_test_case(self.rings, self.atoms, self.basic, one_point=False, two_point=True))

	def _run_test_case(self, molecules, addons=None, linkers=None, one_point=True, two_point=False):
		_molecules = map(lambda s: {'smiles': s}, sample(molecules, min(self.num_repeats, len(molecules))))
		_addons = (
			map(lambda s: {'smiles': s}, sample(addons, min(self.num_repeats, len(addons))))
			if addons else None)
		_linkers = (
			map(lambda s: {'smiles': s}, sample(linkers, min(self.num_repeats, len(linkers))))
			if linkers else None)
		config_fn = create_config(_molecules, _addons, _linkers, one_point, two_point)
		client = CombinationsClient.from_file(config_fn)
		self.files_to_remove.append(config_fn)
		self.paths_to_remove.append(client._config.output.path)
		addons_fn = client._config.config_dict.get('addons', None)
		linkers_fn = client._config.config_dict.get('linkers', None)
		if addons_fn:
			self.files_to_remove.append(addons_fn)
		if linkers_fn:
			self.files_to_remove.append(linkers_fn)
		client.process()
		return self._all_unique(client._config)

	def _all_unique(self, run_config):
		status = True
		path = run_config.output.path

		for fn in map(lambda x: path + os.sep + x, filter(lambda x: x.endswith('smi'), os.listdir(path))):
			handlers = map(MoleculeHandler, IOUtils.read_smi(fn))
			status = status and len(handlers) == len(get_unique_mols(map(lambda x: x.mol.graph, handlers)))

		return status


def create_config(molecules, addons=None, linkers=None, one_point=True, two_point=False):
	temp_fn = tempfile.mkstemp()[1]
	temp_pn = tempfile.mkdtemp()
	molecules_fn = tempfile.mkstemp()[1]
	IOUtils.write_json(molecules_fn, molecules)
	json_config = {
		'molecules': molecules_fn, 
		'output': {'path': temp_pn},
		'attach': {'one_point': one_point, 'two_point': two_point},
		'merge': {'one_point': one_point, 'two_point': two_point}}
	if addons:
		addons_fn = tempfile.mkstemp()[1]
		IOUtils.write_json(addons_fn, addons)
		json_config['addons'] = addons_fn
	if linkers:
		linkers_fn = tempfile.mkstemp()[1]
		IOUtils.write_json(linkers_fn, linkers)
		json_config['linkers'] = linkers_fn
	IOUtils.write_json(temp_fn, json_config)
	return temp_fn


if __name__ == '__main__':
	unittest.main()
