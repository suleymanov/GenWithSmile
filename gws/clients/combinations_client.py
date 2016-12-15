"""
Module for various combinations generation.
"""

import os
import sys
from collections import namedtuple
from json import dumps, loads
from multiprocessing import Pool

from gws.core import Molecule
from gws.core import Modifier

Config = namedtuple('Config', ['X', 'Y', 'L', 'attach', 'merge', 'output', 'numproc'])
Modification = namedtuple('Modification', ['one_point', 'two_point'])


class CombinationsClient(object):
	def __init__(self, config):
		CombinationsClient._validate(config)
		self._first_list = config.X
		self._second_list = config.Y
		self._linkers_list = config.L

		self._attach = config.attach
		self._merge = config.merge
		self._output = config.output
		self._numproc = config.numproc
		if not os.path.exists(self._output['path']):
			os.mkdir(self._output['path'])

	def process(self):
		pack_i = lambda i: (
			i, self._first_list, self._second_list, self._linkers_list,
			self._attach, self._merge, self._output)

		for i_start in xrange(0, len(self._first_list), self._numproc):
			pool = Pool(processes=self._numproc)
			pool.map(_process, map(pack_i, xrange(i_start, i_start + self._numproc)))
			pool.close()
			pool.join()

	@staticmethod
	def from_file(fn):
		"""
		Initialize client from file.
		:param fn: str
		:return: CombinationsClient
		"""
		return CombinationsClient(_read_config(fn))

	@staticmethod
	def _validate(config):
		# TODO: implement smarter exceptions handling
		# TODO: implement incompatible cases handling
		assert config.X is not None
		assert config.Y is not None
		assert config.L is not None
		assert all(len(item['indexes']) == 2 for item in config.L)
		assert config.attach is not None
		assert config.merge is not None
		assert config.output is not None
		assert config.numproc > 0
		assert len(config.X) > 0 and len(config.Y) > 0
		assert (
			config.attach.one_point or config.attach.two_point or
			config.merge.one_point or config.merge.two_point)
		if config.L:
			assert config.attach.one_point or config.merge.one_point
		assert 'path' in config.output
		assert 'file_name_alias' in config.output
		assert 'max_entries_per_file' in config.output


def _update_results(results, i, j, output):
	max_entries = output['max_entries_per_file']
	for ind in xrange(0, len(results) / max_entries + 1):
		fn = '{}{}{}_{}_{}_pack_{}.smi'.format(
			output['path'], os.sep, output['file_name_alias'], i, j, ind + 1)
		with open(fn, 'w') as f:
			f.write('\n'.join(map(
				lambda mol: mol.get_smiles(),
				results[ind*max_entries:(ind+1)*max_entries])) + '\n')


def _process(params):
	i = params[0]
	first_list = params[1]
	second_list = params[2]
	linkers_list = params[3]
	attach_settings = params[4]
	merge_settings = params[5]
	output = params[6]
	if i >= len(first_list):
		return

	hosts = [first_list[i]]
	if linkers_list:
		new_hosts = []
		while len(hosts) > 0:
			modifier = Modifier(hosts.pop())
			for linker_item in linkers_list:
				linker_mol = Molecule.from_smiles(linker_item['linker'])
				positions = linker_item['indexes']
				linker_hosts = []
				if merge_settings.one_point:
					linker_hosts.extend(modifier.merge_explicit(linker_mol, [positions[0]]))
				if attach_settings.one_point:
					linker_hosts.extend(modifier.attach_explicit(linker_mol, [positions[0]]))
				for host in linker_hosts:
					host.reset_next_positions()
					next_pos = host.get_new_index(positions[1])
					host.add_next_positions([next_pos])
				new_hosts.extend(linker_hosts)
		hosts = new_hosts[::]
		for host in hosts:
			host.update_positions()
	for j in xrange(len(second_list)):
		mol_2 = second_list[j]
		results = []
		for host in hosts:
			modifier = Modifier(host)
			results.extend(
				modifier.attach(mol_2, attach_settings.one_point, attach_settings.two_point))
			results.extend(
				modifier.merge(mol_2, merge_settings.one_point, merge_settings.two_point))
		if len(results) > 0:
			_update_results(results, i, j, output)


def _read_config(fn):
	assert os.path.exists(fn)
	with open(fn) as f:
		config_data = loads(f.read())
	return _init_config(config_data)


def _read_smiles(fn):
	assert os.path.exists(fn)
	with open(fn) as f:
		return f.read().splitlines()


def _read_json(fn):
	assert os.path.exists(fn)
	with open(fn) as f:
		return loads(f.read())


def _read_molecules(entry, allowed_atoms):
	fmt = entry['format']
	file_name = entry['file_name']
	if fmt == 'smiles':
		return map(
			lambda x: Molecule.from_smiles(x, allowed_symbols=allowed_atoms),
			_read_smiles(file_name))
	elif fmt == 'json':
		return map(
			lambda x: Molecule(mol_smiles=x['smiles'], positions=x.get('positions', None)),
			_read_json(file_name))
	else:
		raise ValueError('Unknown format.')


def _init_config(config_data):
	linkers_file_name = config_data.get('linkers', None)
	allowed_atoms = config_data.get('allowed_atoms', [])
	return Config(
		X=_read_molecules(config_data.get('first_list', ''), allowed_atoms),
		Y=_read_molecules(config_data.get('second_list', ''), allowed_atoms),
		L=_read_json(linkers_file_name) if linkers_file_name else [],
		attach=Modification(
			one_point=config_data['attach'].get('one_point', True),
			two_point=config_data['attach'].get('two_point', False)),
		merge=Modification(
			one_point=config_data['merge'].get('one_point', True),
			two_point=config_data['merge'].get('two_point', False)),
		output=config_data.get('output_options', {}),
		numproc=config_data.get('numproc', 1))
