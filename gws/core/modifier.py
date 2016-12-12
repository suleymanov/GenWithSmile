"""
Module for modifying molecule.
"""

from algs import get_unique_coords
from modifier_factory import ModifierFactory


class Modifier(object):
	def __init__(self, molecule):
		"""
		Define molecule and coordinates set for modification.
		"""
		self.mol = molecule
		self._one_point_coords = Modifier._define_coords_set(
			self.mol.positions,
			get_unique_coords(self.mol, num_points=1))
		self._two_point_coords = Modifier._define_coords_set(
			self.mol.positions,
			get_unique_coords(self.mol, num_points=2))

	def attach(self, addon_mol, one_point, two_point):
		"""
		Generate all possible attaches.
		:param addon_mol: molecule.Molecule
		:param one_point: bool
		:param two_point: bool
		:return: list of molecule.Molecule
		"""
		results = []
		if one_point:
			results.extend(self._combine(addon_mol, ModifierFactory.attach, 1))
		if two_point:
			results.extend(self._combine(addon_mol, ModifierFactory.attach, 2))
		return results

	def merge(self, addon_mol, one_point, two_point):
		"""
		Generate all possible merges.
		:param addon_mol: molecule.Molecule
		:param one_point: bool
		:param two_point: bool
		:return: list of molecule.Molecule
		"""
		results = []
		if one_point:
			results.extend(self._combine(addon_mol, ModifierFactory.merge, 1))
		if two_point:
			results.extend(self._combine(addon_mol, ModifierFactory.merge, 2))
		return results

	def attach_explicit(self, addon_mol, addon_coord):
		if not 1 <= len(addon_coord) <= 2:
			raise ValueError('Three and more point modifications await their implementation.')
		results = []
		for coord in self._get_coords_set(num_points=len(addon_coord)):
			results.extend(ModifierFactory.attach(self.mol, addon_mol, coord, addon_coord))
		return results

	def merge_explicit(self, addon_mol, addon_coord):
		if not 1 <= len(addon_coord) <= 2:
			raise ValueError('Three and more point modifications await their implementation.')
		results = []
		for coord in self._get_coords_set(num_points=len(addon_coord)):
			results.extend(ModifierFactory.merge(self.mol, addon_mol, coord, addon_coord))
		return results

	def _combine(self, addon_mol, mod_func, num_points):
		results = []
		for coord in self._get_coords_set(num_points):
			for addon_coord in Modifier._define_coords_set(
				addon_mol.positions, get_unique_coords(addon_mol, num_points)):
				results.extend(mod_func(self.mol, addon_mol, coord, addon_coord))
		return results

	def _get_coords_set(self, num_points):
		if num_points == 1:
			return self._one_point_coords
		if num_points == 2:
			return self._two_point_coords
		# TODO: smarter error handling
		raise ValueError('Three and more point modifications await their implementation.')

	@staticmethod
	def _define_coords_set(atom_ids, unique_coords, num_points=1):
		# TODO: make smarter implementation and errors handling
		assert num_points > 0

		if num_points > 2:
			raise ValueError('Three and more point modifications await their implementation.')
		return (
			filter(lambda coord: coord[0] in atom_ids, unique_coords)
			if num_points == 1 else
			filter(lambda coord: coord[0] in atom_ids and coord[1] in atom_ids, unique_coords))
