# Tests for gws.src.core.modifier module

import unittest

from rdkit import Chem

from gws.src.core.algs import rdkitmol2graph
from gws.src.core.modifier import Attacher, Merger


class CoreModifierAttacherTests(unittest.TestCase):
	def setUp(self):
		self.modifier = Attacher

	def test_ctor_1(self):
		"""
		Check that class ctor succeeds/fails when appropriate (benzene as input).
		Case 1: should fail if there are no inputs to the ctor.
		Case 2-5: should not fail; check that:
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		with self.assertRaises(StandardError):
			self.modifier()

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 0)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 2, 4])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

	def test_ctor_2(self):
		"""
		Check that class ctor succeeds/fails when appropriate (piperazine as input).
		Case 1-4: should not fail; check that:
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 4)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 2)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 2, 4])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

	def test_ctor_3(self):
		"""
		Check that class ctor succeeds/fails when appropriate (chloromethanediamine as input).
		Case 1-4: should not fail; check that: 
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		mol_smiles = 'NC(N)Cl'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 4)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 1)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 2)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 1, 2])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

	def test_correct_buddy_1(self):
		"""
		Check that modification succeeds when applied to the modification of the same type.
		Case 1: attach cyclohexane to the attacher of another cyclohexane.
		Case 2: attach cyclohexane to the attacher of benzene.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		modifier.combine(modifier)

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		other_modifier = self.modifier(mol_rdkit, mol_graph)
		modifier.combine(other_modifier)

	def test_incorrect_buddy_1(self):
		"""
		Check that modification fails when applied to the modification of another type.
		Case 1: attach cyclohexane to the merger of another cyclohexane.
		Case 2: attach benzene to the merger of another benzene.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)

		wrong_modifier = Merger(mol_rdkit, mol_graph)
		with self.assertRaises(StandardError):
			modifier.combine(wrong_modifier)

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		wrong_modifier = Merger(mol_rdkit, mol_graph)
		with self.assertRaises(StandardError):
			modifier.combine(wrong_modifier)

	def test_combine_1(self):
		"""
		Perform cyclohexane attach to another cyclohexane. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 2)
		self.assertTrue(len(modifier.one_point_coords) == 1)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 4)
		self.assertTrue(len(modifier.two_point_coords) == 1)

	def test_combine_2(self):
		"""
		Perform benzene attach to another benzene. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 1)
		self.assertTrue(len(modifier.one_point_coords) == 1)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 1)
		self.assertTrue(len(modifier.two_point_coords) == 1)

	def test_combine_3(self):
		"""
		Perform piperazine attach to another piperazine. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 5)
		self.assertTrue(len(modifier.one_point_coords) == 2)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 10)
		self.assertTrue(len(modifier.two_point_coords) == 2)

	def test_combine_4(self):
		"""
		Perform piperazine attach to cyclohexane. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		other_modifier = self.modifier(mol_rdkit, mol_graph)
		other_modifier._init_coords(one_point=True)
		one_p_res = modifier.combine(other_modifier)
		self.assertTrue(len(one_p_res) == 3)
		self.assertTrue(len(modifier.one_point_coords) == 2)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		self.assertTrue(len(other_modifier.one_point_coords) == 1)
		self.assertTrue(len(other_modifier.two_point_coords) == 0)

		other_modifier._init_coords(two_point=True)
		two_p_res = modifier.combine(other_modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 6)
		self.assertTrue(len(modifier.two_point_coords) == 2)
		self.assertTrue(len(other_modifier.two_point_coords) == 1)


class CoreModifierMergerTests(unittest.TestCase):
	def setUp(self):
		self.modifier = Merger

	def test_ctor_1(self):
		"""
		Check that class ctor succeeds/fails when appropriate (benzene as input).
		Case 1: should fail if there are no inputs to the ctor.
		Case 2-5: should not fail; check that:
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		with self.assertRaises(StandardError):
			self.modifier()

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 0)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 2, 4])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

	def test_ctor_2(self):
		"""
		Check that class ctor succeeds/fails when appropriate (piperazine as input).
		Case 1-4: should not fail; check that:
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 6)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 4)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 2)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		
		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 2, 4])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		
	def test_ctor_3(self):
		"""
		Check that class ctor succeeds/fails when appropriate (chloromethanediamine as input).
		Case 1-4: should not fail; check that:
			- available positions were initialized appropriately
			- no 1/2-point orbits were initialized yet.
		"""
		mol_smiles = 'NC(N)Cl'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		self.assertTrue(len(modifier.pos) == 4)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		
		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['C'])
		self.assertTrue(len(modifier.pos) == 1)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		
		modifier = self.modifier(mol_rdkit, mol_graph, atoms=['N'])
		self.assertTrue(len(modifier.pos) == 2)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		
		modifier = self.modifier(mol_rdkit, mol_graph, positions=[0, 1, 2])
		self.assertTrue(len(modifier.pos) == 3)
		self.assertTrue(len(modifier.next_pos) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 0)
		self.assertTrue(len(modifier.two_point_coords) == 0)

	def test_correct_buddy_1(self):
		"""
		Check that modification succeeds when applied to the modification of the same type.
		Case 1: merge cyclohexane to the merger of another cyclohexane.
		Case 2: merge cyclohexane to the merger of benzene.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		modifier.combine(modifier)

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		other_modifier = self.modifier(mol_rdkit, mol_graph)
		modifier.combine(other_modifier)

	def test_incorrect_buddy_1(self):
		"""
		Check that modification fails when applied to the modification of another type.
		Case 1: merge cyclohexane to the attacher of another cyclohexane.
		Case 2: merge benzene to the attacher of another benzene.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)

		wrong_modifier = Attacher(mol_rdkit, mol_graph)
		with self.assertRaises(StandardError):
			modifier.combine(wrong_modifier)

		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		wrong_modifier = Attacher(mol_rdkit, mol_graph)
		with self.assertRaises(StandardError):
			modifier.combine(wrong_modifier)

	def test_combine_1(self):
		"""
		Perform cyclohexane merge to another cyclohexane. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 1)
		self.assertTrue(len(modifier.one_point_coords) == 1)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 2)
		self.assertTrue(len(modifier.two_point_coords) == 1)

	def test_combine_2(self):
		"""
		Perform benzene merge to another benzene. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'c1ccccc1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 0)
		self.assertTrue(len(modifier.one_point_coords) == 1)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 3)
		self.assertTrue(len(modifier.two_point_coords) == 1)

	def test_combine_3(self):
		"""
		Perform piperazine merge to another piperazine. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(modifier)
		self.assertTrue(len(one_p_res) == 2)
		self.assertTrue(len(modifier.one_point_coords) == 2)
		self.assertTrue(len(modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 6)
		self.assertTrue(len(modifier.two_point_coords) == 2)

	def test_combine_4(self):
		"""
		Perform piperazine merge to cyclohexane. Check that:
			- overall number of results is correct
			- 1/2-point orbits are initialized correctly.
		"""
		mol_smiles = 'N1CCNCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		modifier = self.modifier(mol_rdkit, mol_graph)
		mol_smiles = 'C1CCCCC1'
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		other_modifier = self.modifier(mol_rdkit, mol_graph)
		one_p_res = modifier.combine(other_modifier)
		self.assertTrue(len(one_p_res) == 2)
		self.assertTrue(len(modifier.one_point_coords) == 2)
		self.assertTrue(len(modifier.two_point_coords) == 0)
		self.assertTrue(len(other_modifier.one_point_coords) == 1)
		self.assertTrue(len(other_modifier.two_point_coords) == 0)

		two_p_res = modifier.combine(other_modifier, one_point=False, two_point=True)
		self.assertTrue(len(two_p_res) == 3)
		self.assertTrue(len(modifier.two_point_coords) == 2)
		self.assertTrue(len(other_modifier.two_point_coords) == 1)


if __name__ == '__main__':
	unittest.main()
