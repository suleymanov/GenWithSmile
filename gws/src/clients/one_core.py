"""
Module for one core molecule modification.
"""

import os
from multiprocessing import Pool
from datetime import datetime

from gws.src.io import read_one_core_config as read_config
from gws.src.io import write_config
from gws.src.core import MoleculeHandler
from gws.src.core import get_unique_mols


class OneCoreClient(object):
	def __init__(self, config):
		self._config = config
		self._handlers = [MoleculeHandler(
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
			self._perform_iteration(it)
			self._update_results(it_num)

	def _perform_iteration(self, iteration):
		iter_handlers = self._handlers[::]
		for _ in xrange(iteration.max):
			handlers = iter_handlers[::]
			results = []
			pack_i = lambda i: (i, handlers, iteration)

			for i_start in xrange(0, len(handlers), self._config.numthreads):
				pool = Pool(processes=self._config.numthreads)
				pool_results = pool.map(
					_process, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
				pool.close()
				pool.join()
				results.extend(reduce(
					lambda res, item: res + item, filter(lambda x: x, pool_results), []))
			
			iter_handlers = results[::]
			self._iter_results.extend(results)

		self._iter_results = OneCoreClient._filter_non_unique(self._iter_results)

		self._handlers = self._iter_results[::]
		for handler in self._handlers:
			handler.update_modifiers_positions()

	def _update_results(self, it_num):
		max_entries = self._config.output.max
		for i in xrange(len(self._iter_results) / max_entries + 1):
			fn = '{}{}{}_{}_iter_{}_pack_{}.smi'.format(
				self._config.output.path, os.sep, self._config.output.alias, 
				self._dt_str, it_num + 1, i + 1)
			with open(fn, 'w') as f:
				f.write('\n'.join(map(
					lambda handler: handler.mol.smiles,
					self._iter_results[i*max_entries:(i+1)*max_entries])) + '\n')
		self._iter_results = []

	@staticmethod
	def from_file(fn):
		return OneCoreClient(read_config(fn))

	@staticmethod
	def _filter_non_unique(handlers):
		return map(lambda ind: handlers[ind], get_unique_mols(map(lambda x: x.mol.graph, handlers)))


def _process(params):
	i = params[0]
	handlers = params[1]
	iteration = params[2]

	if i >= len(handlers):
		return

	handler = handlers[i]
	return (
		reduce(
			lambda res, x: res + handler.attach(
				MoleculeHandler(
					mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
				iteration.attach.one_point, iteration.attach.two_point),
			iteration.attach.addons, []) + 
		reduce(
			lambda res, x: res + handler.merge(
				MoleculeHandler(
					mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
				iteration.merge.one_point, iteration.merge.two_point),
			iteration.merge.addons, []))
