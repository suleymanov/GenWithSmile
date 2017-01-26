"""
Module for combinations generation between two sets of fragments.
(With linkers if provided.)
"""

import os
from multiprocessing import Pool
from datetime import datetime

from gws.src.io import read_combinations_config as read_config
from gws.src.io import write_config
from gws.src.core import MoleculeHandler
from gws.src.core import get_unique_mols


class CombinationsClient(object):
	def __init__(self, config):
		self._config = config
		self._mol_handlers = map(
			lambda x: MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.molecules)
		self._addon_handlers = map(
			lambda x: MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.addons)
		self._linker_handlers = map(
			lambda x: MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.linkers)
		for i, handler in enumerate(self._mol_handlers):
			handler.create_mol_graph('handler_{}'.format(i + 1))
		for i, handler in enumerate(self._addon_handlers):
			handler.create_mol_graph('{}_{}'.format(
				'handler' if self._config.samelist else 'addon', i + 1))
		for i, handler in enumerate(self._linker_handlers):
			handler.create_mol_graph('linker_{}'.format(i + 1))
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		self.write_config(config_fn)

	def process(self):
		pack_i = lambda i: (
			i, self._mol_handlers, self._addon_handlers, self._linker_handlers, self._config,
			self._dt_str)

		for i_start in xrange(0, len(self._mol_handlers), self._config.numthreads):
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
	mol_handlers = params[1]
	addon_handlers = params[2]
	linker_handlers = params[3]
	config = params[4]
	dt_str = params[5]
	shift = i if config.samelist else 0
	addon_handlers = addon_handlers[shift:]

	if i >= len(mol_handlers):
		return

	handlers = [mol_handlers[i]]
	if linker_handlers:
		handler = handlers.pop()
		
		for linker_handler in linker_handlers:
			handlers.extend(
				handler.attach(linker_handler, config.attach.one_point, config.attach.two_point) +
				handler.merge(linker_handler, config.merge.one_point, config.merge.two_point))
	for j in xrange(len(addon_handlers)):
		addon_handler = addon_handlers[j]
		results = []
		for handler in handlers:
			results.extend(
				handler.attach(addon_handler, config.attach.one_point, config.attach.two_point) +
				handler.merge(addon_handler, config.merge.one_point, config.merge.two_point))
			if len(results) > 0:
				_update_results(_filter_non_unique(results), i, j + shift, config.output, dt_str)


def _filter_non_unique(handlers):
	return map(lambda ind: handlers[ind], get_unique_mols(map(lambda x: x.mol_graph, handlers)))


def _update_results(results, i, j, output, dt_str):
	max_entries = output.max
	for ind in xrange(0, len(results) / max_entries + 1):
		fn = '{}{}{}_{}_{}_{}_pack_{}.smi'.format(
			output.path, os.sep, output.alias, dt_str, i, j, ind + 1)
		with open(fn, 'w') as f:
			f.write('\n'.join(map(
				lambda x: x.mol_smiles,
				results[ind*max_entries:(ind+1)*max_entries])) + '\n')
