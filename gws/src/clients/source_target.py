"""
Module for linkers generation between source and target molecules.
"""

import os
from multiprocessing import Pool
from datetime import datetime

from gws.src.io import read_source_target_config as read_config
from gws.src.io import write_config
from gws.src.core import Modifier, get_unique_mols


class SourceTargetClient(object):
	def __init__(self, config):
		self._config = config
		self._modifiers = [Modifier(
			mol_smiles=self._config.source.smiles, atoms=self._config.source.atoms,
			attach_pos=self._config.source.attach_pos, merge_pos=self._config.source.merge_pos)]
		self._target_modifier = Modifier(
			mol_smiles=self._config.target.smiles, atoms=self._config.target.atoms,
			attach_pos=self._config.target.attach_pos, merge_pos=self._config.target.merge_pos)
		self._iter_results = []
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		write_config(config_fn, self._config)

	def process(self):
		for it_num, it in enumerate(self._config.iterations):
			self._perform_iteration(it)
			self._update_results(it_num)

	def _perform_iteration(self, iteration):
		iter_modifiers = self._modifiers[::]
		src_w_linkers = []
		pack_i = lambda i: (i, iter_modifiers, iteration)

		for i_start in xrange(0, len(iter_modifiers), self._config.numthreads):
			pool = Pool(processes=self._config.numthreads)
			pool_results = pool.map(
				_add_linkers, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
			pool.close()
			pool.join()
			src_w_linkers.extend(reduce(
				lambda res, item: res + item, filter(lambda x: x, pool_results), []))

		for modifier in src_w_linkers:
			modifier.update_positions()

		pack_i = lambda i: (i, src_w_linkers, iteration, self._target_modifier)
		self._modifiers = []

		for i_start in xrange(0, len(src_w_linkers), self._config.numthreads):
			pool = Pool(processes=self._config.numthreads)
			pool_results = pool.map(
				_extend, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
			pool.close()
			pool.join()
			self._iter_results.extend(reduce(
				lambda res, item: res + item, filter(lambda x: x, pool_results), []))
			for i in xrange(i_start, i_start + self._config.numthreads):
				if pool_results[i - i_start]:
					self._modifiers.append(src_w_linkers[i])

		self._iter_results = get_unique_mols(self._iter_results)
		self._modifiers = get_unique_mols(self._modifiers)

	def _update_results(self, it_num):
		max_entries = self._config.output.max
		for i in xrange(len(self._iter_results) / max_entries + 1):
			fn = '{}{}{}_{}_iter_{}_pack_{}.smi'.format(
				self._config.output.path, os.sep, self._config.output.alias, 
				self._dt_str, it_num + 1, i + 1)
			with open(fn, 'w') as f:
				f.write('\n'.join(map(
					lambda modifier: modifier.get_smiles(),
					self._iter_results[i*max_entries:(i+1)*max_entries])) + '\n')
		self._iter_results = []

	@staticmethod
	def from_file(fn):
		return SourceTargetClient(read_config(fn))


def _add_linkers(params):
	i = params[0]
	modifiers = params[1]
	iteration = params[2]

	if i >= len(modifiers):
		return

	modifier = modifiers[i]
	results = []
	for addon in iteration.attach.addons:
		addon_modifier = Modifier(
			mol_smiles=addon.smiles, atoms=addon.atoms,
			attach_pos=addon.attach_pos, merge_pos=addon.merge_pos)
		results.extend(modifier.attach(
			addon_modifier, iteration.attach.one_point, iteration.attach.two_point))
	for addon in iteration.merge.addons:
		addon_modifier = Modifier(
			mol_smiles=addon.smiles, atoms=addon.atoms,
			attach_pos=addon.attach_pos, merge_pos=addon.merge_pos)
		results.extend(modifier.merge(
			addon_modifier, iteration.merge.one_point, iteration.merge.two_point))
	return results


def _extend(params):
	i = params[0]
	modifiers = params[1]
	iteration = params[2]
	target_modifier = params[3]

	if i >= len(modifiers):
		return

	modifier = modifiers[i]
	return (
		modifier.attach(target_modifier, iteration.attach.one_point, iteration.attach.two_point) +
		modifier.merge(target_modifier, iteration.merge.one_point, iteration.merge.two_point))
