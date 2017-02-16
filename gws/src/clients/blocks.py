"""
Module for iterative generating from ordered blocks.
"""

import os
from datetime import datetime

from gws.src.io import read_blocks_config as read_config
from gws.src.io import write_config
from gws.src.core import MoleculeHandler
from gws.src.core import get_unique_mols, conn_index


class BlocksClient(object):
	def __init__(self, config):
		self._config = config
		self._addon_handlers = map(
			lambda x: MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.blocks)
		self._handlers = map(lambda x: [x], self._addon_handlers)
		self._results = [[]] * len(self._handlers)
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(
			dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		self._write_config(config_fn)

	def process(self):
		for _ in xrange(self._config.max):
			for i_h in xrange(len(self._handlers)):
				results = []
				while len(self._handlers[i_h]) > 0:
					curr_handler = self._handlers[i_h].pop()
					for addon_handler in self._addon_handlers[:(i_h + 1)]:
						results.extend(curr_handler.attach(addon_handler))
				for res in results:
					res.update_modifiers_positions()
				self._results[i_h] = self._results[i_h][::] + results[::]
				self._handlers[i_h] = results[::]
		self._update_results()

	def _update_results(self):
		from rdkit import Chem

		for i, res in enumerate(self._results):
			fn = '{}{}{}_{}_{}.smi'.format(
				self._config.output.path, os.sep, self._config.output.alias, self._dt_str, i + 1)
			filtered = filter(
				lambda handler: (all(map(
					lambda s: handler.mol.rdkit.HasSubstructMatch(Chem.MolFromSmarts(s)),
					self._config.patterns.include))
				and all(map(
					lambda s: not handler.mol.rdkit.HasSubstructMatch(Chem.MolFromSmarts(s)),
					self._config.patterns.exclude))),
				BlocksClient._filter_non_unique(res))
			filtered = filter(
				lambda handler: all(map(
					lambda bound: bound.low<=conn_index(handler.mol.rdkit, bound.expon)<=bound.high, 
					self._config.profile)),
				filtered)
			with open(fn, 'w') as f:
				f.write('\n'.join(map(lambda handler: handler.mol.smiles, filtered)))

	def _write_config(self, fn):
		config_dict = {
			'molecules': self._config.blocks,
			'depth': self._config.max,
			'output': self._config.output,
			'numthreads': self._config.numthreads
		}
		if 'patterns' in self._config.config_dict:
			config_dict['patterns'] = self._config.config_dict['patterns']
		write_config(fn, config_dict)

	@staticmethod
	def from_file(fn):
		return BlocksClient(read_config(fn))

	@staticmethod
	def _filter_non_unique(handlers):
		return map(lambda ind: handlers[ind], get_unique_mols(map(lambda x: x.mol.graph, handlers)))
