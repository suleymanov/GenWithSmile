import os
from collections import namedtuple

from simplejson import dumps

from gws.src.io.utils import IOUtils
from validation_factory import ValidationFactory
from reader_tools import ReaderTools


OneCoreSettings = namedtuple(
	'OneCoreSettings', 
	['core', 'iterations', 'output', 'patterns', 'numthreads', 'config_dict'])
SourceTargetSettings = namedtuple(
	'SourceTargetSettings', 
	['source', 'target', 'iterations', 'output', 'patterns', 'numthreads', 'config_dict'])
CombinationsSettings = namedtuple(
	'CombinationsSettings',
	['molecules', 'addons', 'linkers', 'attach', 'merge', 'output', 'patterns', 'numthreads', 
	'samelist', 'config_dict'])
BlocksSettings = namedtuple(
	'BlocksSettings', ['blocks', 'max', 'output', 'patterns', 'profile', 'numthreads', 'config_dict'])
SynthpathSettings = namedtuple('SynthpathSettings', ['blocks', 'target', 'output', 'rules'])


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
			patterns=ReaderTools.read_patterns_settings(config_dict.get('patterns', {}), IOUtils),
			numthreads=config_dict.get('numthreads', 1),
			config_dict=config_dict
		)

	@staticmethod
	def _validate_config(config):
		ValidationFactory.validate_molecule(config.core)
		ValidationFactory.validate_iterations(config.iterations)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_threading(config.numthreads)
		ValidationFactory.validate_patterns(config.patterns.include)
		ValidationFactory.validate_patterns(config.patterns.exclude)


class BlocksConfigReader(object):
	def __init__(self, config_dict):
		BlocksConfigReader._check_config_dict(config_dict)
		config = BlocksConfigReader._create_config(config_dict)
		BlocksConfigReader._validate_config(config)
		self.config = config
		if not os.path.exists(self.config.output.path):
			os.mkdir(self.config.output.path)

	@staticmethod
	def from_file(fn):
		return BlocksConfigReader(IOUtils.read_json(fn))

	@staticmethod
	def _check_config_dict(config_dict):
		assert 'molecules' in config_dict, 'No building blocks settings found.'

	@staticmethod
	def _create_config(config_dict):
		return BlocksSettings(
			blocks=map(ReaderTools.read_molecule, config_dict['molecules']),
			max=config_dict.get('depth', 1),
			output=ReaderTools.read_output(config_dict.get('output', {})),
			patterns=ReaderTools.read_patterns_settings(config_dict.get('patterns', {}), IOUtils),
			profile=map(ReaderTools.read_bound_settings, config_dict.get('profile', [])),
			numthreads=config_dict.get('numthreads', 1),
			config_dict=config_dict
		)

	@staticmethod
	def _validate_config(config):
		map(ValidationFactory.validate_molecule, config.blocks)
		map(ValidationFactory.validate_bound, config.profile)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_patterns(config.patterns.include)
		ValidationFactory.validate_patterns(config.patterns.exclude)
		ValidationFactory.validate_threading(config.numthreads)
		exps = map(lambda x: x.expon, config.profile)
		assert len(exps) == len(set(exps)), 'All exponents values must be unique.'
		assert config.max > 0, 'Should be at least one iteration.'


class SynthpathConfigReader(object):
	def __init__(self, config_dict):
		SynthpathConfigReader._check_config_dict(config_dict)
		config = SynthpathConfigReader._create_config(config_dict)
		SynthpathConfigReader._validate_config(config)
		self.config = config
		if not os.path.exists(self.config.output.path):
			os.mkdir(self.config.output.path)

	@staticmethod
	def from_file(fn):
		return SynthpathConfigReader(IOUtils.read_json(fn))

	@staticmethod
	def _check_config_dict(config_dict):
		assert 'blocks' in config_dict, 'No building blocks settings found.'
		assert 'target' in config_dict, 'No synthesis target found.'

	@staticmethod
	def _create_config(config_dict):
		return SynthpathSettings(
			blocks=map(ReaderTools.read_molecule, config_dict.get('blocks', [])),
			target=ReaderTools.read_molecule(config_dict['target']),
			output=ReaderTools.read_output(config_dict.get('output', {})),
			rules=[
				'*C(=O)N(*)',    	# amide
				'*C(=O)O*',         # ester
				'*N(*)(*)',         # amine
				'*N(*)C(=O)N(*)*',  # urea
				'*O*',              # ether
				'*C(*)=C(*)*',      # olefin
				'*[N+](*)(*)*',     # quaternary nitrogen
				'nC',               # aromatic nitrogen-aliphatic carbon
				'[N;R][C;!R]',      # lactam nitrogen-aliphatic carbon
				'c-c',              # aromatic carbon-aromatic carbon
				'*N(*)S(=O)(=O)*'   # sulphonamide
			])

	@staticmethod
	def _validate_config(config):
		map(ValidationFactory.validate_molecule, config.blocks)
		ValidationFactory.validate_molecule, config.target
		ValidationFactory.validate_output(config.output)


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
			patterns=ReaderTools.read_patterns_settings(config_dict.get('patterns', {}), IOUtils),
			numthreads=config_dict.get('numthreads', 1),
			config_dict=config_dict
		)

	@staticmethod
	def _validate_config(config):
		ValidationFactory.validate_molecule(config.source)
		ValidationFactory.validate_molecule(config.target)
		ValidationFactory.validate_output(config.output)
		ValidationFactory.validate_threading(config.numthreads)
		ValidationFactory.validate_iterations(config.iterations)
		ValidationFactory.validate_patterns(config.patterns.include)
		ValidationFactory.validate_patterns(config.patterns.exclude)
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
			patterns=ReaderTools.read_patterns_settings(config_dict.get('patterns', {}), IOUtils),
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
		ValidationFactory.validate_patterns(config.patterns.include)
		ValidationFactory.validate_patterns(config.patterns.exclude)


def read_one_core_config(fn):
	reader = OneCoreConfigReader.from_file(fn)
	return reader.config


def read_source_target_config(fn):
	reader = SourceTargetConfigReader.from_file(fn)
	return reader.config


def read_combinations_config(fn):
	reader = CombinationsConfigReader.from_file(fn)
	return reader.config


def read_blocks_config(fn):
	reader = BlocksConfigReader.from_file(fn)
	return reader.config


def read_synthpath_config(fn):
	reader = SynthpathConfigReader.from_file(fn)
	return reader.config
