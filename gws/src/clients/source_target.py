"""
Module for linkers generation between source and target molecules.
"""

import os
from multiprocessing import Pool
from datetime import datetime

from gws.src.io import read_source_target_config as read_config
from gws.src.io import write_config
from gws.src.core import MoleculeHandler
from gws.src.core import get_unique_mols


class SourceTargetClient(object):
	def __init__(self, config):
		self._config = config
		self._handlers = [MoleculeHandler(
			mol_smiles=self._config.source.smiles, atoms=self._config.source.atoms,
			attach_pos=self._config.source.attach_pos, merge_pos=self._config.source.merge_pos)]
		self._target_handler = MoleculeHandler(
			mol_smiles=self._config.target.smiles, atoms=self._config.target.atoms,
			attach_pos=self._config.target.attach_pos, merge_pos=self._config.target.merge_pos)
		self._iter_results = []
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		self._write_config(config_fn)

	def process(self):
		for it_num, it in enumerate(self._config.iterations):
			self._perform_iteration(it)
			self._update_results(it_num)

	def _perform_iteration(self, iteration):
		iter_handlers = self._handlers[::]
		src_w_linkers = []
		pack_i = lambda i: (i, iter_handlers, iteration)

		for i_start in xrange(0, len(iter_handlers), self._config.numthreads):
			pool = Pool(processes=self._config.numthreads)
			pool_results = pool.map(
				_add_linkers, map(pack_i, xrange(i_start, i_start + self._config.numthreads)))
			pool.close()
			pool.join()
			src_w_linkers.extend(reduce(
				lambda res, item: res + item, filter(lambda x: x, pool_results), []))

		for handler in src_w_linkers:
			handler.update_modifiers_positions()

		pack_i = lambda i: (i, src_w_linkers, iteration, self._target_handler)
		self._handlers = []

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
					self._handlers.append(src_w_linkers[i])

		self._iter_results = SourceTargetClient._filter_non_unique(self._iter_results)
		self._handlers = SourceTargetClient._filter_non_unique(self._handlers)

	def _update_results(self, it_num):
		from rdkit import Chem

		max_entries = self._config.output.max
		for i in xrange(len(self._iter_results) / max_entries + 1):
			fn = '{}{}{}_{}_iter_{}_pack_{}.smi'.format(
				self._config.output.path, os.sep, self._config.output.alias, 
				self._dt_str, it_num + 1, i + 1)
			with open(fn, 'w') as f:
				f.write('\n'.join(map(
					lambda handler: handler.mol.smiles,
					filter(
						lambda handler: (all(map(
							lambda s: handler.mol.rdkit.HasSubstructMatch(Chem.MolFromSmarts(s)),
							self._config.patterns.include))
						and all(map(
							lambda s: not handler.mol.rdkit.HasSubstructMatch(Chem.MolFromSmarts(s)),
							self._config.patterns.exclude))),
						self._iter_results[i*max_entries:(i+1)*max_entries]))) + '\n')
		self._iter_results = []

	def _write_config(self, fn):
		config_dict = {
			'source': self._config.source,
			'target': self._config.target,
			'iterations': self._config.iterations,
			'output': self._config.output,
			'numthreads': self._config.numthreads
		}
		if 'patterns' in self._config.config_dict:
			config_dict['patterns'] = self._config.config_dict['patterns']
		write_config(fn, config_dict)

	@staticmethod
	def from_file(fn):
		return SourceTargetClient(read_config(fn))

	@staticmethod
	def _filter_non_unique(handlers):
		return map(lambda ind: handlers[ind], get_unique_mols(map(lambda x: x.mol.graph, handlers)))


def _add_linkers(params):
	i = params[0]
	handlers = params[1]
	iteration = params[2]

	if i >= len(handlers):
		return

	handler = handlers[i]
	return (
		reduce(lambda res, x: res + handler.attach(
			MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			iteration.attach.one_point, iteration.attach.two_point),
			iteration.attach.addons, []) +
		reduce(lambda res, x: res + handler.merge(
			MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			iteration.attach.one_point, iteration.attach.two_point),
			iteration.merge.addons, []))


def _extend(params):
	i = params[0]
	handlers = params[1]
	iteration = params[2]
	target_handler = params[3]

	if i >= len(handlers):
		return

	handler = handlers[i]
	return (
		handler.attach(target_handler, iteration.attach.one_point, iteration.attach.two_point) +
		handler.merge(target_handler, iteration.merge.one_point, iteration.merge.two_point))
