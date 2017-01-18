"""
Module for combinations generation between two sets of fragments.
(With linkers if provided.)
"""

import os
from multiprocessing import Pool
from datetime import datetime

from gws.src.io import read_combinations_config as read_config
from gws.src.io import write_config
from gws.src.core import Modifier, get_unique_mols


class CombinationsClient(object):
	def __init__(self, config):
		self._config = config
		self._mol_modifiers = map(
			lambda x: Modifier(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.molecules)
		self._addon_modifiers = map(
			lambda x: Modifier(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.addons)
		self._linker_modifiers = map(
			lambda x: Modifier(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.linkers)
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		self.write_config(config_fn)

	def process(self):
		pack_i = lambda i: (
			i, self._mol_modifiers, self._addon_modifiers, self._linker_modifiers, self._config,
			self._dt_str)

		for i_start in xrange(0, len(self._mol_modifiers), self._config.numthreads):
			pool = Pool(processes=self._config.numthreads)
			pool.map(_process, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
			pool.close()
			pool.join()

	@staticmethod
	def from_file(fn):
		return CombinationsClient(read_config(fn))

	def write_config(self, fn):
		config_dict = {
			'molecules': self._config.config_dict['molecules'],
			'attach': self._config.attach,
			'merge': self._config.merge,
			'output': self._config.output,
			'numthreads': self._config.numthreads,
		}
		addons_fn = self._config.config_dict.get('addons', None)
		linkers_fn = self._config.config_dict.get('linkers', None)
		if addons_fn:
			config_dict['addons'] = addons_fn
		if linkers_fn:
			config_dict['linkers'] = linkers_fn
		write_config(fn, config_dict)


def _process(params):
	i = params[0]
	mol_modifiers = params[1]
	addon_modifiers = params[2]
	linker_modifiers = params[3]
	config = params[4]
	dt_str = params[5]
	shift = i if config.samelist else 0
	addon_modifiers = addon_modifiers[shift:]

	if i >= len(mol_modifiers):
		return

	modifiers = [mol_modifiers[i]]
	if linker_modifiers:
		modifier = modifiers.pop()
		for linker_modifier in linker_modifiers:
			modifiers.extend(
				modifier.attach(linker_modifier, config.attach.one_point, config.attach.two_point) +
				modifier.merge(linker_modifier, config.merge.one_point, config.merge.two_point))
		for modifier in modifiers:
			modifier.update_positions()
	for j in xrange(len(addon_modifiers)):
		addon_modifier = addon_modifiers[j]
		results = []
		for modifier in modifiers:
			results.extend(
				modifier.attach(addon_modifier, config.attach.one_point, config.attach.two_point) +
				modifier.merge(addon_modifier, config.merge.one_point, config.merge.two_point))
			if len(results) > 0:
				_update_results(get_unique_mols(results), i, j + shift, config.output, dt_str)


def _update_results(results, i, j, output, dt_str):
	max_entries = output.max
	for ind in xrange(0, len(results) / max_entries + 1):
		fn = '{}{}{}_{}_{}_{}_pack_{}.smi'.format(
			output.path, os.sep, output.alias, dt_str, i, j, ind + 1)
		with open(fn, 'w') as f:
			f.write('\n'.join(map(
				lambda mol: mol.get_smiles(),
				results[ind*max_entries:(ind+1)*max_entries])) + '\n')
