# Base class for arbitrary modification type

from abc import ABCMeta, abstractmethod
from collections import namedtuple

from rdkit import Chem
from networkx import Graph, MultiGraph

from gws.src.core.algs import get_unique_coords


CombinationResult = namedtuple('CombinationResult', ['result_mol_rdkit', 'pos', 'next_pos', 'inds_map', 'ind'])


bonds_def = {
	1.0: Chem.rdchem.BondType.SINGLE,
	2.0: Chem.rdchem.BondType.DOUBLE,
	3.0: Chem.rdchem.BondType.TRIPLE,
	4.0: Chem.rdchem.BondType.QUADRUPLE,
	5.0: Chem.rdchem.BondType.QUINTUPLE
}


class BaseModifier(object):
	__metaclass__ = ABCMeta

	def __init__(self, mol_rdkit, positions=None, atoms=None):
		"""
		New version of an abstract Modifier instance (without explicitly defining molecular graph).
		:param mol_rdkit: Chem.rdchem.Mol
		:param positions: list of int
		:param atoms: list of str
		"""
		self._mol_rdkit = mol_rdkit
		self._mol_graph = None
		all_pos = (
			map(
				lambda x: x.GetIdx(),
				filter(lambda x: x.GetSymbol() in atoms, self._mol_rdkit.GetAtoms()))
			if atoms else list(xrange(self._mol_rdkit.GetNumAtoms())))
		self.pos = filter(lambda x: x in positions, all_pos) if positions else all_pos
		self.next_pos = []
		self.one_point_coords = []
		self.two_point_coords = []

	def set_mol_graph(self, mol_graph):
		self._mol_graph = mol_graph

	def combine(self, other_modifier, one_point=True, two_point=False):
		"""
		Generate all possible combinations of two fragments according to selected modifier.
		:param other_modifier: BaseModifier
		:param one_point: bool
		:param two_point: bool
		:return: list of CombinationResult
		"""
		assert type(self) == type(other_modifier)

		results = []
		if one_point:
			for ind, coord in enumerate(self.get_coords(one_point=one_point)):
				for other_coord in other_modifier.get_coords(one_point=one_point):
					results.extend(self._combine_one_point(other_modifier, coord, other_coord, ind))
		if two_point:
			for ind, coord in enumerate(self.get_coords(two_point=two_point)):
				for other_coord in other_modifier.get_coords(two_point=two_point):
					results.extend(self._combine_two_points(other_modifier, coord, other_coord, ind))
		return results

	def get_coords(self, one_point=False, two_point=False):
		"""
		Return 1 or 2 point orbits (initialize when queried for the first time).
		:param one_point: bool
		:param two_point: bool
		:return: list of list of int
		"""
		if one_point:
			if not self.one_point_coords:
				self._init_coords(one_point=one_point)
			return self.one_point_coords
		if two_point:
			if not self.two_point_coords:
				self._init_coords(two_point=two_point)
			return self.two_point_coords
		raise KeyError('Only 1/2-point coordinates are available.')

	def add_next_pos(self, next_pos):
		"""
		Add positions for future modifications.
		:param next_pos: list of int
		"""
		self.next_pos = list(set(self.next_pos + next_pos))

	def update_pos(self):
		"""
		Update positions for future modifications.
		"""
		self.pos = self.next_pos[::]
		self.next_pos = []

	def _init_coords(self, one_point=False, two_point=False):
		if one_point:
			self.one_point_coords = get_unique_coords(
				self._mol_graph, 1, map(lambda x: [x], self.pos))
		if two_point:
			self.two_point_coords = filter(
				lambda x: x[0] in self.pos and x[1] in self.pos,
				map(lambda x: (x.GetBeginAtomIdx(), x.GetEndAtomIdx()), self._mol_rdkit.GetBonds()))

	def _pack_result(self, other_modifier, result_mol, inds_map, coord, other_coord, ind):
		curr_pos = filter(lambda x: x not in coord, self.pos)
		next_pos = map(
			lambda x: inds_map[x], 
			filter(lambda x: x not in other_coord and x in other_modifier.pos, inds_map))
		next_pos += self.next_pos
		return CombinationResult(result_mol, curr_pos, next_pos, inds_map, ind)

	@abstractmethod
	def _combine_one_point(self, other_modifier, coord, other_coord, ind):
		pass

	@abstractmethod
	def _combine_two_points(self, other_modifier, coord, other_coord, ind):
		pass
