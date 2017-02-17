"""
Module for best synthesis path construction.
"""

import os
from datetime import datetime
from heapq import heappush, heappop
from json import dumps

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

from gws.src.io import read_synthpath_config as read_config
from gws.src.io import write_config
from gws.src.core import MoleculeHandler


class SynthpathClient(object):
	def __init__(self, config):
		self._config = config
		self._blocks = map(
			lambda x: MoleculeHandler(
				mol_smiles=x.smiles, atoms=x.atoms, attach_pos=x.attach_pos, merge_pos=x.merge_pos),
			self._config.blocks)
		self._target = MoleculeHandler(
			mol_smiles=self._config.target.smiles, atoms=self._config.target.atoms,
			attach_pos=self._config.target.attach_pos, merge_pos=self._config.target.merge_pos)
		self._target_fp = FingerprintMols.FingerprintMol(self._target.mol.rdkit)
		self._rules = map(Chem.MolFromSmarts, self._config.rules)
		self._pq = []
		for block_handler in self._blocks:
			fp = FingerprintMols.FingerprintMol(block_handler.mol.rdkit)
			diff = 1.0 - DataStructs.FingerprintSimilarity(fp, self._target_fp)
			heappush(self._pq, (0 + diff, {
				'handler': block_handler, 'path': [], 'init': block_handler, 'addons': [], 'steps': 0}))
		self._best_path = None
		dt = datetime.now()
		self._dt_str = '{}_{}_{}_{}_{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
		config_fn = '{}{}{}_{}_config.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		write_config(config_fn, self._config)

	def process(self):
		results = []
		while len(self._pq) > 0:
			handler_data = heappop(self._pq)[1]
			handler = handler_data['handler']
			path = handler_data['path']
			init = handler_data['init']
			for i, block_handler in enumerate(self._blocks):
				block_results = handler.attach(block_handler)
				for j, res in enumerate(block_results):
					rule = self._get_rule(handler.mol.rdkit, res.mol.rdkit)
					if rule:
						if res.mol.rdkit.GetNumAtoms() > self._target.mol.rdkit.GetNumAtoms():
							continue
						fp = FingerprintMols.FingerprintMol(res.mol.rdkit)
						init_fp = FingerprintMols.FingerprintMol(init.mol.rdkit)
						diff_0 = 1.0 - DataStructs.FingerprintSimilarity(fp, init_fp)
						diff = 1.0 - DataStructs.FingerprintSimilarity(fp, self._target_fp)
						res_path = path + [rule]
						addons = handler_data['addons'] + [Chem.MolToSmiles(block_handler.mol.rdkit)]
						steps = handler_data['steps'] + 1
						res.update_modifiers_positions()
						new_item = {
							'handler': res, 'path': res_path, 'init': init, 'addons': addons, 'steps': steps}
						if diff == 0.0:
							results.append(new_item)
						else:
							heappush(self._pq, (diff_0 + diff, new_item))
		if len(results) > 0:
			best_res = sorted(results, key=lambda x: x['steps'])[0]
			self._best_path = {
				'initial': Chem.MolToSmiles(best_res['init'].mol.rdkit),
				'result': Chem.MolToSmiles(best_res['handler'].mol.rdkit),
				'synthpath': map(
					lambda x: {'pattern': x[0], 'addon': x[1]}, 
					zip(best_res['path'], best_res['addons']))
			}
		else:
			self._best_path = {'msg': 'No paths found.'}
		self._update_results()

	def _get_rule(self, first_mol, result_mol):
		rules_matches = []
		for rule in self._rules:
			match = filter(
				lambda m: (
					any(map(lambda x: x < first_mol.GetNumAtoms(), m)) and 
					any(map(lambda x: x >= first_mol.GetNumAtoms(), m))),
				result_mol.GetSubstructMatches(rule))
			if len(match) > 0:
				rules_matches.append(Chem.MolToSmarts(rule))
		if len(rules_matches) > 0:
			return rules_matches[0]

	def _update_results(self):
		fn = '{}{}{}_{}.json'.format(
			self._config.output.path, os.sep, self._config.output.alias, self._dt_str)
		with open(fn, 'w') as f:
			f.write(dumps(self._best_path, indent=4))

	@staticmethod
	def from_file(fn):
		return SynthpathClient(read_config(fn))
