import os
from collections import namedtuple

from simplejson import dumps

from gws.src.io.utils import IOUtils
from validation_factory import ValidationFactory
from reader_tools import ReaderTools


OneCoreSettings = namedtuple(
	'OneCoreSettings', 
	['core', 'iterations', 'output', 'numthreads'])
SourceTargetSettings = namedtuple(
	'SourceTargetSettings', 
	['source', 'target', 'iterations', 'output', 'numthreads'])
CombinationsSettings = namedtuple(
	'CombinationsSettings',
	['molecules', 'addons', 'linkers', 'attach', 'merge', 'output', 'numthreads', 'samelist', 'config_dict'])


class OneCoreConfigReader(object):
	def __init__(self, config_dict):
		OneCoreConfigReader._check_config_dict(config_dict)
		config = OneCoreConfigReader._create_config(config_dict)
		OneCoreConfigReader._validate_config(config)
		self.config = config
		if not os.path.exists(self.config.output.path):
			os.mkdir(self.config.output.path)

	@staticmethod
	def from_file(fn):
		return OneCoreConfigReader(IOUtils.read_json(fn))

	@staticmethod
	def _check_config_dict(config_dict):
		assert 'core' in config_dict, 'No molecule settings found.'
		assert 'iterations' in config_dict, 'No iterations settings found.'

	@staticmethod
	def _create_config(config_dict):
		return OneCoreSettings(
			core=ReaderTools.read_molecule(config_dict['core']),
			iterations=map(ReaderTools.read_iteration, config_dict['iterations']),
			output=ReaderTools.read_output(config_dict.get('output', {})),
			numthreads=config_dict.get('numthreads', 1)
		)

	@staticmethod
	def _validate_config(config):
		ValidationFactory.validate_molecule(config.core)
		ValidationFactory.validate_iterations(config.iterations)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_threading(config.numthreads)


class SourceTargetConfigReader(object):
	def __init__(self, config_dict):
		SourceTargetConfigReader._check_config_dict(config_dict)
		config = SourceTargetConfigReader._create_config(config_dict)
		SourceTargetConfigReader._validate_config(config)
		self.config = config
		if not os.path.exists(self.config.output.path):
			os.mkdir(self.config.output.path)

	@staticmethod
	def from_file(fn):
		return SourceTargetConfigReader(IOUtils.read_json(fn))

	@staticmethod
	def _check_config_dict(config_dict):
		assert 'source' in config_dict, 'No source molecule settings found.'
		assert 'target' in config_dict, 'No target molecule settings found.'
		assert 'iterations' in config_dict, 'No iterations settings found.'

	@staticmethod
	def _create_config(config_dict):
		return SourceTargetSettings(
			source=ReaderTools.read_molecule(config_dict['source']),
			target=ReaderTools.read_molecule(config_dict['target']),
			iterations=map(ReaderTools.read_iteration, config_dict['iterations']),
			output=ReaderTools.read_output(config_dict.get('output', {})),
			numthreads=config_dict.get('numthreads', 1)
		)

	@staticmethod
	def _validate_config(config):
		ValidationFactory.validate_molecule(config.source)
		ValidationFactory.validate_molecule(config.target)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_threading(config.numthreads)
		ValidationFactory.validate_iterations(config.iterations)
		assert all(it.max == 1 for it in config.iterations), \
			'Handling more than one addon awaits for implementation.'


class CombinationsConfigReader(object):
	def __init__(self, config_dict):
		CombinationsConfigReader._check_config_dict(config_dict)
		config = CombinationsConfigReader._create_config(config_dict)
		CombinationsConfigReader._validate_config(config)
		self.config = config
		if not os.path.exists(self.config.output.path):
			os.mkdir(self.config.output.path)

	@staticmethod
	def from_file(fn):
		return CombinationsConfigReader(IOUtils.read_json(fn))

	@staticmethod
	def _check_config_dict(config_dict):
		assert 'molecules' in config_dict, 'No molecules file name found.'

	@staticmethod
	def _create_config(config_dict):
		mols_file_name = config_dict['molecules']
		adds_file_name = config_dict.get('addons', mols_file_name)
		linkers_file_name = config_dict.get('linkers', '')

		return CombinationsSettings(
			molecules=map(ReaderTools.read_molecule, IOUtils.read_json(mols_file_name)),
			addons=map(ReaderTools.read_molecule, IOUtils.read_json(adds_file_name)),
			linkers=(
				map(ReaderTools.read_molecule, IOUtils.read_json(linkers_file_name))
				if linkers_file_name else []),
			attach=ReaderTools.read_modification(config_dict.get('attach', {})),
			merge=ReaderTools.read_modification(config_dict.get('merge', {})),
			output=ReaderTools.read_output(config_dict.get('output', {})),
			numthreads=config_dict.get('numthreads', 1),
			samelist=mols_file_name==adds_file_name,
			config_dict=config_dict
		)

	@staticmethod
	def _validate_config(config):
		map(ValidationFactory.validate_molecule, config.molecules)
		map(ValidationFactory.validate_molecule, config.addons)
		map(ValidationFactory.validate_molecule, config.linkers)
		ValidationFactory.validate_modification(config.attach)
		ValidationFactory.validate_modification(config.merge)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_threading(config.numthreads)


def read_one_core_config(fn):
	reader = OneCoreConfigReader.from_file(fn)
	return reader.config


def read_source_target_config(fn):
	reader = SourceTargetConfigReader.from_file(fn)
	return reader.config


def read_combinations_config(fn):
	reader = CombinationsConfigReader.from_file(fn)
	return reader.config


def write_config(fn, config):
	with open(fn, 'w') as f:
		f.write(dumps(config, indent=4))