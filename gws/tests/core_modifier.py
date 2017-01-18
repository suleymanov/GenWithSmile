import unittest
from collections import namedtuple

from rdkit import Chem

from gws.src.core.algs import (
	is_isomorph,
	get_unique_coords,
	get_unique_mols,
	rdkitmol2graph,
	get_nonisomorphic_positions
)
from gws.src.core.modifier import Modifier

FakeMolecule = namedtuple('FakeMolecule', ['rdkit', 'graph'])
FakeModifier = namedtuple('FakeModifier', ['mol'])


class CoreModifierTests(unittest.TestCase):
	def test_1_modifier_ctor(self):
		with self.assertRaises(StandardError):
			Modifier()
		
		Modifier(mol_smiles='c1ccccc1')
		
		modifier = Modifier(mol_rdkit=Chem.MolFromSmiles('c1ccccc1'))
		self.assertTrue(len(modifier._attach_pos) == 6)
		self.assertTrue(len(modifier._merge_pos) == 6)
		self.assertTrue(len(modifier._next_attach_pos) == 0)
		self.assertTrue(len(modifier._next_merge_pos) == 0)
		
		modifier = Modifier(mol_rdkit=Chem.MolFromSmiles('c1ccccc1'), atoms=['N'])
		self.assertTrue(len(modifier._attach_pos) == 0)
		self.assertTrue(len(modifier._merge_pos) == 0)

		modifier = Modifier(mol_smiles='c1ccccc1', attach_pos=[0, 2, 4], merge_pos=[1, 3, 5])
		self.assertTrue(len(modifier._attach_pos) == 3)
		self.assertTrue(len(modifier._merge_pos) == 3)
		self.assertTrue(all(x not in modifier._merge_pos for x in modifier._attach_pos))

	def test_2_modifier_create_coords(self):
		modifier = Modifier(mol_smiles='c1ccccc1')
		modifier.create_attach_coords(True, True)
		modifier.create_merge_coords(True, True)
		self.assertTrue(len(modifier.get_attach_coords(1)) == 1)
		self.assertTrue(len(modifier.get_attach_coords(2)) == 1)
		self.assertTrue(len(modifier.get_merge_coords(1)) == 1)
		self.assertTrue(len(modifier.get_merge_coords(2)) == 1)

		modifier = Modifier(mol_smiles='c1ccncc1')
		modifier.create_attach_coords(True, True)
		modifier.create_merge_coords(True, True)
		self.assertTrue(len(modifier.get_attach_coords(1)) == 4)
		self.assertTrue(len(modifier.get_attach_coords(2)) == 3)
		self.assertTrue(len(modifier.get_merge_coords(1)) == 4)
		self.assertTrue(len(modifier.get_merge_coords(2)) == 3)

		modifier = Modifier(mol_smiles='c1ccncc1', atoms=['C'])
		modifier.create_attach_coords(True, True)
		modifier.create_merge_coords(True, True)
		self.assertTrue(len(modifier.get_attach_coords(1)) == 3)
		self.assertTrue(len(modifier.get_attach_coords(2)) == 2)
		self.assertTrue(len(modifier.get_merge_coords(1)) == 3)
		self.assertTrue(len(modifier.get_merge_coords(2)) == 2)

	def test_3_attach(self):
		host_modifier = Modifier(mol_smiles='c1ccccc1')
		host_modifier.create_attach_coords(True, True)
		self.assertTrue(len(host_modifier.attach(host_modifier, one_point=False, two_point=False)) == 0)
		self.assertTrue(len(host_modifier.attach(host_modifier, one_point=True, two_point=False)) == 1)
		self.assertTrue(len(host_modifier.attach(host_modifier, one_point=False, two_point=True)) == 1)
		self.assertTrue(len(host_modifier.attach(host_modifier, one_point=True, two_point=True)) == 2)

	def test_4_merge(self):
		host_modifier = Modifier(mol_smiles='C1CCCCC1')
		host_modifier.create_merge_coords(True, True)
		self.assertTrue(len(host_modifier.merge(host_modifier, one_point=False, two_point=False)) == 0)
		self.assertTrue(len(host_modifier.merge(host_modifier, one_point=True, two_point=False)) == 1)
		self.assertTrue(len(host_modifier.merge(host_modifier, one_point=False, two_point=True)) == 2)
		self.assertTrue(len(host_modifier.merge(host_modifier, one_point=True, two_point=True)) == 3)


if __name__ == '__main__':
	unittest.main()
