"""
Module for modifying one molecule.
"""

from collections import namedtuple

from rdkit import Chem

from algs import get_unique_coords
from algs import rdkitmol2graph
from modifier_factory import ModifierFactory


Molecule = namedtuple('Molecule', ['smiles', 'rdkit', 'graph'])


class Modifier(object):
	def __init__(self, mol_smiles=None, mol_rdkit=None, atoms=None, attach_pos=None, merge_pos=None):
		"""
		Initialize molecule modifier with appropriate attach/merge positions.
		:param mol_smiles: str
		:param mol_rdkit: Chem.rdchem.Mol
		:param atoms: list of str
		:param attach_pos: list of int
		:param merge_pos: list of int
		:return: None
		"""
		self._init_molecule(mol_smiles, mol_rdkit)
		self._init_positions(atoms, attach_pos, merge_pos)
		self._one_point_attach_coords = []
		self._two_point_attach_coords = []
		self._one_point_merge_coords = []
		self._two_point_merge_coords = []

	def create_attach_coords(self, one_point, two_point):
		"""
		Create coordinates that will be used in attach modifications.
		:param one_point: bool
		:param two_point: bool
		:return: None
		"""
		self._one_point_attach_coords = Modifier._define_coords_set(
			self._attach_pos,
			get_unique_coords(self.mol, 1, self._attach_pos)) if one_point else []
		self._two_point_attach_coords = Modifier._define_coords_set(
			self._attach_pos,
			get_unique_coords(self.mol, 2, self._attach_pos)) if two_point else []

	def create_merge_coords(self, one_point, two_point):
		"""
		Create coordinates that will be used in merge modifications.
		:param one_point: bool
		:param two_point: bool
		:return: None
		"""
		self._one_point_merge_coords = Modifier._define_coords_set(
			self._merge_pos,
			get_unique_coords(self.mol, 1, self._merge_pos)) if one_point else []
		self._two_point_merge_coords = Modifier._define_coords_set(
			self._merge_pos,
			get_unique_coords(self.mol, 2, self._merge_pos)) if two_point else []

	def attach(self, addon_modifier, one_point, two_point):
		"""
		Perform 1 or 2 point attach of an addon modifier.
		:param addon_modifier: Modifier
		:param one_point: bool
		:param two_point: bool
		:return: list of Modifier
		"""
		results = []
		if one_point:
			results.extend(self._combine(addon_modifier, 1, attach=True))
		if two_point:
			results.extend(self._combine(addon_modifier, 2, attach=True))
		return results

	def merge(self, addon_modifier, one_point, two_point):
		"""
		Perform 1 or 2 point merge of an addon modifier.
		:param addon_modifier: Modifier
		:param one_point: bool
		:param two_point: bool
		:return: list of Modifier
		"""
		results = []
		if one_point:
			results.extend(self._combine(addon_modifier, 1, merge=True))
		if two_point:
			results.extend(self._combine(addon_modifier, 2, merge=True))
		return results

	def get_attach_coords(self, num_points):
		if num_points == 1:
			if not self._one_point_attach_coords:
				self.create_attach_coords(one_point=True, two_point=False)
			return self._one_point_attach_coords
		if num_points == 2:
			if not self._two_point_attach_coords:
				self.create_attach_coords(one_point=False, two_point=True)
			return self._two_point_attach_coords
		raise KeyError('Only 1/2-point coordinates are available.')

	def get_merge_coords(self, num_points):
		if num_points == 1:
			if not self._one_point_merge_coords:
				self.create_merge_coords(one_point=True, two_point=False)
			return self._one_point_merge_coords
		if num_points == 2:
			if not self._two_point_merge_coords:
				self.create_merge_coords(one_point=False, two_point=True)
			return self._two_point_merge_coords
		raise KeyError('Only 1/2-point coordinates are available.')

	def get_smiles(self):
		return self.mol.smiles

	def add_next_positions(self, next_attach_pos, next_merge_pos):
		"""
		Add positions that will be used in the further steps.
		:param next_attach_pos: list of int
		:param next_merge_pos: list of int
		:return: None
		"""
		self._next_attach_pos = list(set(self._next_attach_pos + next_attach_pos))
		self._next_merge_pos = list(set(self._next_merge_pos + next_merge_pos))

	def update_positions(self):
		"""
		Update positions to prepare modifier to the further steps.
		:return: None
		"""
		self._attach_pos = self._next_attach_pos[::]
		self._merge_pos = self._next_merge_pos[::]
		self._next_attach_pos = []
		self._next_merge_pos = []

	def create_new_modifier(
		self, result_mol, addon_modifier, addon_inds_map, rm_attach_pos=None, rm_merge_pos=None):
		"""
		Create new modifier from modification results.
		:param result_mol: Chem.rdchem.Mol
		:param addon_inds_map: dict  # {addon_old_index: result_new_index}
		:param rm_attach_pos: list of int
		:param rm_merge_pos: list of int
		:return: Modifier
		"""
		attach_pos = self._attach_pos[::]
		if rm_attach_pos is not None:
			attach_pos = filter(lambda x: x not in rm_attach_pos, attach_pos)
		merge_pos = self._merge_pos[::]
		if rm_merge_pos is not None:
			merge_pos = filter(lambda x: x not in rm_merge_pos, merge_pos)
		modifier = Modifier(mol_rdkit=result_mol, attach_pos=attach_pos, merge_pos=merge_pos)
		modifier.add_next_positions(
			next_attach_pos=map(
				lambda k: addon_inds_map[k], 
				filter(lambda k: k in addon_modifier._attach_pos, addon_inds_map)),
			next_merge_pos=map(
				lambda k: addon_inds_map[k], 
				filter(lambda k: k in addon_modifier._merge_pos, addon_inds_map)))
		return modifier

	def _init_molecule(self, mol_smiles, mol_rdkit):
		# TODO: implement initialization from networkx graph
		if mol_smiles is not None:
			_mol_smiles = mol_smiles
			_mol_rdkit = Chem.MolFromSmiles(_mol_smiles)
		elif mol_rdkit is not None:
			_mol_rdkit = mol_rdkit
			_mol_smiles = Chem.MolToSmiles(_mol_rdkit)
		else:
			raise ValueError('Should initialize either from SMILES or rdkit graph.')
		_mol_graph = rdkitmol2graph(_mol_rdkit)
		self.mol = Molecule(_mol_smiles, _mol_rdkit, _mol_graph)

	def _init_positions(self, atoms=None, attach_pos=None, merge_pos=None):
		all_pos = (
			map(
				lambda at: at.GetIdx(),
				filter(lambda at: at.GetSymbol() in atoms, self.mol.rdkit.GetAtoms()))
			if atoms else list(xrange(self.mol.rdkit.GetNumAtoms())))
		self._attach_pos = filter(lambda x: x in all_pos, attach_pos) if attach_pos else all_pos
		self._merge_pos = filter(lambda x: x in all_pos, merge_pos) if merge_pos else all_pos
		self._next_attach_pos = []
		self._next_merge_pos = []

	def _combine(self, addon_modifier, num_points, attach=False, merge=False):
		results = []
		if attach:
			for coord in self.get_attach_coords(num_points):
				for addon_coord in addon_modifier.get_attach_coords(num_points):
					results.extend(ModifierFactory.attach(self, addon_modifier, coord, addon_coord))
		if merge:
			for coord in self.get_merge_coords(num_points):
				for addon_coord in addon_modifier.get_merge_coords(num_points):
					results.extend(ModifierFactory.merge(self, addon_modifier, coord, addon_coord))
		return results

	@staticmethod
	def _define_coords_set(atom_inds, unique_coords, num_points=1):
		assert 0 < num_points < 3, 'Only 1/2-point modifications are implemented.'
		return (
			filter(lambda coord: coord[0] in atom_inds, unique_coords)
			if num_points == 1 else
			filter(lambda coord: coord[0] in atom_inds and coord[1] in atom_inds, unique_coords))
