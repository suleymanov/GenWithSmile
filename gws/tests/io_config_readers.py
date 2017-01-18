import unittest
import tempfile

from gws.src.io import IOUtils
from gws.src.io import read_one_core_config
from gws.src.io import read_source_target_config
from gws.src.io import read_combinations_config

from utils import clean_files


class OneCoreConfigReaderTests(unittest.TestCase):
	def setUp(self):
		self.fn_config = tempfile.mkstemp()[1]
		self.read_config = read_one_core_config

	def tearDown(self):
		clean_files([self.fn_config])

	def test_minimal_valid(self):
		config_dict = {
			'core': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)  # should not fail

	def test_full_valid(self):
		config_dict = {
			'core': {
				'smiles': 'C1CC1', 
				'atoms': ['C'], 
				'attach_pos': [0, 1, 2], 
				'merge_pos': [0, 1, 2]
			},
			'iterations': [
				{
					'attach': {
						'addons': [{
							'smiles': 'C', 
							'attach_pos': [0], 
							'merge_pos': [], 
							'atoms': ['C']
						}],
						'one_point': True,
						'two_point': False
					},
					'merge': {
						'addons': [{
							'smiles': 'O', 
							'attach_pos': [0], 
							'merge_pos': [], 
							'atoms': ['O']
						}],
						'one_point': True,
						'two_point': False
					},
					'max': 1
				}
			],
			'output': {
				'path': 'one_core_results',
				'alias': 'one_core_result',
				'max': 1000
			},
			'numthreads': 1
		}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)

	def test_minimal_invalid(self):
		config_dict = {}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_core(self):
		config_dict = {'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {'core': {}, 'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_core(self):
		config_dict = {
			'core': {'smiles': 123},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'abc'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_atoms(self):
		config_dict = {
			'core': {'smiles': 'C1CC1', 'atoms': {}},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'C1CC1', 'atoms': [1]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_positions(self):
		config_dict = {
			'core': {'smiles': 'C1CC1', 'attach_pos': {}},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'C1CC1', 'attach_pos': [-1]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'C1CC1', 'attach_pos': ['0']},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'C1CC1', 'attach_pos': [3]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'core': {'smiles': 'C1CC1', 'attach_pos': [0, 0]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_iterations(self):
		config_dict = {'core': {'smiles': 'C1CC1'}}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {'core': {'smiles': 'C1CC1'}, 'iterations': []}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_modifications(self):
		config_dict = {'core': {'smiles': 'C1CC1'}, 'iterations': [{'max': 1}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_addons(self):
		config_dict = {'core': {'smiles': 'C1CC1'}, 'iterations': [{'attach': {'addons': []}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_threading(self):
		config_dict = {
			'core': {'smiles': 'C1CC1'}, 
			'iterations': [{'attach': {'addons': [{'smiles': 'C1CC1'}]}}], 
			'numthreads': 0}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

		import multiprocessing as mp

		config_dict = {
			'core': {'smiles': 'C1CC1'}, 
			'iterations': [{'attach': {'addons': [{'smiles': 'C1CC1'}]}}], 
			'numthreads': mp.cpu_count() + 1}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)


class SourceTargetConfigReaderTests(unittest.TestCase):
	def setUp(self):
		self.fn_config = tempfile.mkstemp()[1]
		self.read_config = read_source_target_config

	def tearDown(self):
		clean_files([self.fn_config])

	def test_minimal_valid(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)

	def test_full_valid(self):
		config_dict = {
			'source': {
				'smiles': 'C1CC1', 
				'atoms': ['C'], 
				'attach_pos': [0, 1, 2], 
				'merge_pos': [0, 1, 2]
			},
			'target': {
				'smiles': 'C1CC1', 
				'atoms': ['C'], 
				'attach_pos': [0, 1, 2], 
				'merge_pos': [0, 1, 2]
			},
			'iterations': [
				{
					'attach': {
						'addons': [{
							'smiles': 'C', 
							'attach_pos': [0], 
							'merge_pos': [], 
							'atoms': ['C']
						}],
						'one_point': True,
						'two_point': False
					},
					'merge': {
						'addons': [{
							'smiles': 'O', 
							'attach_pos': [0], 
							'merge_pos': [], 
							'atoms': ['O']
						}],
						'one_point': True,
						'two_point': False
					},
					'max': 1
				}
			],
			'output': {
				'path': 'source_target_results',
				'alias': 'source_target_result',
				'max': 1000
			},
			'numthreads': 1
		}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)

	def test_minimal_invalid(self):
		config_dict = {}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_source(self):
		config_dict = {
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_target(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_source(self):
		config_dict = {
			'source': {'smiles': 123},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'abc'},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_target(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 123},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'abc'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_atoms(self):
		config_dict = {
			'source': {'smiles': 'C1CC1', 'atoms': {}},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1', 'atoms': [1]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_positions(self):
		config_dict = {
			'source': {'smiles': 'C1CC1', 'attach_pos': {}},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1', 'attach_pos': [-1]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1', 'attach_pos': ['0']},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1', 'attach_pos': [3]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1', 'attach_pos': [0, 0]},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_iterations(self):
		config_dict = {'source': {'smiles': 'C1CC1'}, 'target': {'smiles': 'C1CC1'}}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		config_dict = {'source': {'smiles': 'C1CC1'}, 'target': {'smiles': 'C1CC1'}, 'iterations': []}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_modifications(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'}, 
			'target': {'smiles': 'C1CC1'}, 
			'iterations': [{'max': 1}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_no_addons(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'}, 
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': []}}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_max(self):
		config_dict = {
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C'}]}, 'max': 2}]}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_threading(self):
		config_dict = { 
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C1CC1'}]}}], 
			'numthreads': 0}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

		import multiprocessing as mp

		config_dict = { 
			'source': {'smiles': 'C1CC1'},
			'target': {'smiles': 'C1CC1'},
			'iterations': [{'attach': {'addons': [{'smiles': 'C1CC1'}]}}], 
			'numthreads': mp.cpu_count() + 1}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)


class CombinationsConfigReaderTests(unittest.TestCase):
	def setUp(self):
		self.fn_config = tempfile.mkstemp()[1]
		self.fn_molecules = tempfile.mkstemp()[1]
		self.fn_addons = tempfile.mkstemp()[1]
		self.fn_linkers = tempfile.mkstemp()[1]
		self.read_config = read_combinations_config

	def tearDown(self):
		clean_files([self.fn_config, self.fn_molecules, self.fn_addons, self.fn_linkers])

	def test_minimal_valid(self):
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'atoms': ['C']}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)

	def test_full_valid(self):
		molecules = [{
			'smiles': 'C1CC1', 'atoms': ['C', 'N'], 'attach_pos': [0, 1, 2], 'merge_pos': [0, 1, 2]}]
		addons = [{
			'smiles': 'C1NC1', 'atoms': ['C', 'N'], 'attach_pos': [0, 1, 2], 'merge_pos': [0, 1, 2]}]
		linkers = [{
			'smiles': 'CN(C)C', 'atoms': ['C', 'N'], 'attach_pos': [0, 2, 3], 'merge_pos': [0, 2, 3]}]
		IOUtils.write_json(self.fn_molecules, molecules)
		IOUtils.write_json(self.fn_addons, addons)
		IOUtils.write_json(self.fn_linkers, linkers)
		config_dict = {
			'molecules': self.fn_molecules,
			'addons': self.fn_addons,
			'linkers': self.fn_linkers,
			'attach': {'one_point': True, 'two_point': True},
			'merge': {'one_point': True, 'two_point': True},
			'output': {
				'path': 'combination_results',
				'alias': 'combination_result',
				'max': 1000
			},
			'numthreads': 1
		}
		IOUtils.write_json(self.fn_config, config_dict)
		self.read_config(self.fn_config)

	def test_minimal_invalid(self):
		config_dict = {}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_molecules(self):
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'abc'}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 123}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_atoms(self):
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'atoms': {}}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'atoms': [1]}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_positions(self):
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'attach_pos': {}}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'attach_pos': [-1]}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'attach_pos': ['0']}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'attach_pos': [3]}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1', 'attach_pos': [0, 0]}])
		config_dict = {'molecules': self.fn_molecules}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

	def test_invalid_threading(self):
		IOUtils.write_json(self.fn_molecules, [{'smiles': 'C1CC1'}])
		config_dict = {'molecules': self.fn_molecules, 'numthreads': 0}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)

		import multiprocessing as mp

		config_dict = {'molecules': self.fn_molecules, 'numthreads': mp.cpu_count() + 1}
		IOUtils.write_json(self.fn_config, config_dict)
		with self.assertRaises(AssertionError):
			self.read_config(self.fn_config)


if __name__ == '__main__':
	unittest.main()
