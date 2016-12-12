"""
Module for one core molecule modification.
"""

import os
from collections import namedtuple
from json import dumps, loads

from gws.core import get_unique_mols
from gws.core import Molecule
from gws.core import Modifier

Config = namedtuple('Config', ['host', 'iterations', 'output'])
Iteration = namedtuple('Iteration', ['attach', 'merge', 'max'])
Modification = namedtuple('Modification', ['addons', 'one_point', 'two_point'])


class OneCoreClient(object):
	def __init__(self, config):
		"""
		Initialize client with settings.
		:param config: Config
		:return: None
		"""
		OneCoreClient._validate(config)
		self._hosts = [config.host]
		self._iterations = config.iterations
		self._output = config.output
		self._iter_results = []
		if not os.path.exists(self._output['path']):
			os.mkdir(self._output['path'])

	def process(self):
		"""
		Perform all generation here.
		:return: None
		"""
		for it_num, it in enumerate(self._iterations):
			self._perform_iteration(it)
			self._update_results(it_num)

	def _perform_iteration(self, it):
		iter_hosts = self._hosts[::]
		for _ in xrange(it.max):
			hosts = iter_hosts[::]
			iter_hosts = []
			while len(hosts) > 0:
				modifier = Modifier(hosts.pop())
				for add in it.attach.addons:
					mod_results = modifier.attach(add, it.attach.one_point, it.attach.two_point)
					iter_hosts.extend(mod_results)
				for add in it.merge.addons:
					mod_results = modifier.merge(add, it.merge.one_point, it.merge.two_point)
					iter_hosts.extend(mod_results)
			iter_hosts = get_unique_mols(iter_hosts)
			self._iter_results.extend(iter_hosts)
		
		self._hosts = self._iter_results[::]
		for host in self._hosts:
			host.update_positions()

	def _update_results(self, it_num):
		max_entries = self._output['max_entries_per_file']
		for i in xrange(0, len(self._iter_results) / max_entries + 1):
			fn = '{}{}{}_iter_{}_pack_{}.smi'.format(
				self._output['path'], os.sep, self._output['file_name_alias'], it_num + 1, i + 1)
			with open(fn, 'w') as f:
				f.write('\n'.join(map(
					lambda mol: mol.get_smiles(), 
					self._iter_results[i*max_entries:(i+1)*max_entries])) + '\n')
		self._iter_results = []

	@staticmethod
	def from_file(fn):
		"""
		Initialize client from file.
		:param fn: str
		:return: OneCoreClient
		"""
		return OneCoreClient(_read_config(fn))

	@staticmethod
	def _validate(config):
		# TODO: implement smarter exceptions handling
		# TODO: implement incompatible cases handling
		assert config.host is not None
		assert config.iterations is not None
		assert config.output is not None

		for it in config.iterations:
			assert len(it.attach.addons) + len(it.merge.addons) > 0
			assert (
				it.attach.one_point or it.attach.two_point or 
				it.merge.one_point or it.merge.one_point)
			assert it.max > 0
		assert 'path' in config.output
		assert 'file_name_alias' in config.output
		assert 'max_entries_per_file' in config.output


def _read_config(fn):
	with open(fn) as f:
		config_data = loads(f.read())
	return _init_config(config_data)


def _read_addons(modification, mol_init_func):
	return Modification(
		addons=map(lambda s: mol_init_func(s), modification.get('addons', [])),
		one_point=modification.get('one_point', True),
		two_point=modification.get('two_point', False))


def _init_config(config_data):
	mol_smiles = config_data.get('core', '')
	mol_settings = config_data.get('core_settings', {})
	output = config_data.get('output_options', {})
	iterations = []
	for iteration in config_data.get('iterations', []):
		iterations.append(Iteration(
			attach=_read_addons(iteration['attach'], Molecule.from_smiles),
			merge=_read_addons(iteration['merge'], Molecule.from_smiles),
			max=iteration.get('max', 1)))
	return Config(
		host=Molecule.from_smiles(mol_smiles, mol_settings.get('allowed_atoms', [])),
		iterations=iterations,
		output=output)
