from collections import namedtuple

from rdkit import Chem

from algs import rdkitmol2graph
from modifier import Attacher, Merger


Molecule = namedtuple('Molecule', ['smiles', 'rdkit', 'graph'])


class MoleculeHandler(object):
	def __init__(self, mol_smiles=None, mol_rdkit=None, attach_pos=None, merge_pos=None, atoms=None):
		"""
		Create main molecule handler.
		:param mol_smiles: str
		:param mol_rdkit: Chem.rdchem.Mol
		:param attach_pos: list of int
		:param merge_pos: list of int
		:param atoms: list of str
		"""
		if mol_smiles is not None:
			_mol_smiles = mol_smiles
			_mol_rdkit = Chem.MolFromSmiles(_mol_smiles)
		elif mol_rdkit is not None:
			_mol_rdkit = mol_rdkit
			_mol_smiles = Chem.MolToSmiles(_mol_rdkit)
		else:
			raise ValueError('Can initialize either from SMILES or rdkit representation.')
		_mol_graph = rdkitmol2graph(_mol_rdkit)
		self.mol = Molecule(_mol_smiles, _mol_rdkit, _mol_graph)
		self.attacher = Attacher(self.mol.rdkit, self.mol.graph, attach_pos, atoms)
		self.merger = Merger(self.mol.rdkit, self.mol.graph, merge_pos, atoms)

	def attach(self, addon_handler, one_point=True, two_point=False):
		"""
		Perform all unique attaches of some addon.
		:param addon_handler: MoleculeHandler
		:param one_point: bool
		:param two_point: bool
		:return: list of MoleculeHandler
		"""
		return map(
			lambda x: MoleculeHandler._create_handler(
				x.result_mol_rdkit,
				x.pos,
				self.merger.pos,
				x.next_pos,
				self.merger.next_pos + map(
					lambda k: x.inds_map[k],
					filter(lambda k: k in addon_handler.merger.pos, x.inds_map))),
			self.attacher.combine(addon_handler.attacher, one_point, two_point))

	def merge(self, addon_handler, one_point=True, two_point=False):
		"""
		Perform all unique merges of some addon.
		:param addon_handler: MoleculeHandler
		:param one_point: bool
		:param two_point: bool
		:return: list of MoleculeHandler
		"""
		return map(
			lambda x: MoleculeHandler._create_handler(
				x.result_mol_rdkit,
				self.attacher.pos,
				x.pos,
				self.attacher.next_pos + map(
					lambda k: x.inds_map[k],
					filter(lambda k: k in addon_handler.attacher.pos, x.inds_map)),
				x.next_pos),
			self.merger.combine(addon_handler.merger, one_point, two_point))

	def update_modifiers_positions(self):
		self.attacher.update_pos()
		self.merger.update_pos()

	@staticmethod
	def _create_handler(mol_rdkit, attach_pos, merge_pos, next_attach_pos, next_merge_pos):
		handler = MoleculeHandler(mol_rdkit=mol_rdkit, attach_pos=attach_pos, merge_pos=merge_pos)
		handler.attacher.add_next_pos(next_attach_pos)
		handler.merger.add_next_pos(next_merge_pos)
		return handler
