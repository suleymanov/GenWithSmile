import os
from collections import namedtuple


MoleculeSettings = namedtuple('MoleculeSettings', ['smiles', 'atoms', 'attach_pos', 'merge_pos'])
ModificationSettings = namedtuple('ModificationSettings', ['addons', 'one_point', 'two_point'])
IterationSettings = namedtuple('IterationSettings', ['attach', 'merge', 'max'])
PatternsSettings = namedtuple('PatternsSettings', ['include', 'exclude'])
OutputSettings = namedtuple('OutputSettings', ['path', 'alias', 'max'])
BoundSettings = namedtuple('BoundSettings', ['alias', 'low', 'high'])


class ReaderTools(object):
	@staticmethod
	def read_molecule(mol_dict):
		return MoleculeSettings(
			smiles=mol_dict.get('smiles', None),
			atoms=mol_dict.get('atoms', []),
			attach_pos=mol_dict.get('attach_pos', []),
			merge_pos=mol_dict.get('merge_pos', [])
		)

	@staticmethod
	def read_modification(mod_dict):
		return ModificationSettings(
			addons=map(ReaderTools.read_molecule, mod_dict.get('addons', [])),
			one_point=mod_dict.get('one_point', True),
			two_point=mod_dict.get('two_point', False)
		)

	@staticmethod
	def read_iteration(iter_dict):
		return IterationSettings(
			attach=ReaderTools.read_modification(iter_dict.get('attach', {})),
			merge=ReaderTools.read_modification(iter_dict.get('merge', {})),
			max=iter_dict.get('max', 1)
		)

	@staticmethod
	def read_patterns_settings(patt_dict, io_utils):
		def process_type(key, patt_dict):
			patts = patt_dict.get(key, {}).get('from_list', [])
			fn = patt_dict.get(key, {}).get('from_df', None)
			if fn:
				patts.extend(io_utils.read_df(fn))
			return patts

		return PatternsSettings(
			include=process_type('include', patt_dict), 
			exclude=process_type('exclude', patt_dict))

	@staticmethod
	def read_bound_settings(bound_dict):
		return BoundSettings(
			alias=bound_dict.get('alias', ''),
			low=bound_dict.get('lo', 0),
			high=bound_dict.get('hi', 0))

	@staticmethod
	def read_output(out_dict):
		return OutputSettings(
			path=out_dict.get('path', "gws_results"),
			alias=out_dict.get('alias', 'result'),
			max=out_dict.get('max', 1000)
		)
