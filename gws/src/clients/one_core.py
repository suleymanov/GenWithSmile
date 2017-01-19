"""
Module for one core molecule modification.
"""

import os
from multiprocessing import Pool
from datetime import datetime
import timeit

from gws.src.io import read_one_core_config as read_config
from gws.src.io import write_config
from gws.src.core import Modifier, get_unique_mols


class OneCoreClient(object):
	def __init__(self, config):
		self._config = config
		self._modifiers = [Modifier(
			mol_smiles=self._config.core.smiles, atoms=self._config.core.atoms,
			attach_pos=self._config.core.attach_pos, merge_pos=self._config.core.merge_pos)]
		self._iter_results = []
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		write_config(config_fn, self._config)

	def process(self):
		for it_num, it in enumerate(self._config.iterations):
			print('Start iteration #{}'.format(it_num + 1))
			start_time = timeit.default_timer()
			self._perform_iteration(it)
			print('Iteration #{} time: {}'.format(it_num + 1, timeit.default_timer() - start_time))
			self._update_results(it_num)

	def _perform_iteration(self, iteration):
		iter_modifiers = self._modifiers[::]
		print('\tStart repeats')
		# for _ in xrange(iteration.max):
		for i in xrange(iteration.max):
			print('\tStart repeat #{}'.format(i + 1))
			modifiers = iter_modifiers[::]
			results = []
			pack_i = lambda i: (i, modifiers, iteration)

			print('\tStart generation: {} initial modifiers.'.format(len(modifiers)))
			start_time = timeit.default_timer()
			for i_start in xrange(0, len(modifiers), self._config.numthreads):
				pool = Pool(processes=self._config.numthreads)
				pool_results = pool.map(
					_process, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
				pool.close()
				pool.join()
				results.extend(reduce(
					lambda res, item: res + item, filter(lambda x: x, pool_results), []))
			print('\tGeneration time: {}.'.format(timeit.default_timer() - start_time))
			print('\t{} structures generated.'.format(len(results)))

			print('\tStarted filtering results.')
			start_time = timeit.default_timer()
			results = get_unique_mols(results)
			print('\tFiltering time: {}.'.format(timeit.default_timer() - start_time))
			print('\t{} unique structures.'.format(len(results)))
			iter_modifiers = results[::]
			self._iter_results.extend(results)

		# self._iter_results = get_unique_mols(self._iter_results)
		self._modifiers = self._iter_results[::]
		for host in self._modifiers:
			host.update_positions()

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
		return OneCoreClient(read_config(fn))


def _process(params):
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
