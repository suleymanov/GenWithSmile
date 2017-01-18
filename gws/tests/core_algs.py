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

FakeMolecule = namedtuple('FakeMolecule', ['rdkit', 'graph'])
FakeModifier = namedtuple('FakeModifier', ['mol'])


class CoreAlgsTests(unittest.TestCase):
	def test_1_isomorphism(self):
		# trivial
		mol = rdkitmol2graph(Chem.MolFromSmiles('c1ccccc1'))
		self.assertTrue(is_isomorph(mol, mol))
		mol = rdkitmol2graph(Chem.MolFromSmiles('C1CNCCC1'))
		self.assertTrue(is_isomorph(mol, mol))
		mol = rdkitmol2graph(Chem.MolFromSmiles('C'))
		self.assertTrue(is_isomorph(mol, mol))
		mol = rdkitmol2graph(Chem.MolFromSmiles('CC(N)C'))
		self.assertTrue(is_isomorph(mol, mol))

		# non-trivial
		self.assertTrue(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['c1ccccc1', 'C1=CC=CC=C1']))))
		self.assertTrue(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['C1CNCCC1', 'C1NCCCC1']))))
		self.assertTrue(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['CC(N)C', 'CC(C)N']))))

		# non isomorphs
		self.assertFalse(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['C1=CC=CC=C1', 'C1CCCCC1']))))
		self.assertFalse(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['CC(N)C', 'CN(C)C']))))
		self.assertFalse(is_isomorph(
			*map(rdkitmol2graph, map(Chem.MolFromSmiles, ['C', 'N']))))

	def test_2_unique_coords(self):
		# one trivial ring
		mol_rdkit = Chem.MolFromSmiles('C1CCCCC1')
		mol_graph = rdkitmol2graph(mol_rdkit)
		fake_mol = FakeMolecule(mol_rdkit, mol_graph)
		self.assertTrue(len(get_unique_coords(fake_mol)) == 1)
		self.assertTrue(len(get_unique_coords(fake_mol, num_points=2)) == 1)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(0, 5))) == 1)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(1, 6))) == 1)
		self.assertTrue(get_unique_coords(fake_mol, positions_to_check=range(0, 5))[0][0] == 0)
		self.assertTrue(get_unique_coords(fake_mol, positions_to_check=range(1, 6))[0][0] == 1)
		with self.assertRaises(StandardError):
			get_unique_coords(fake_mol, num_points=3)

		# one ring
		mol_rdkit = Chem.MolFromSmiles('C1CNCCC1')
		mol_graph = rdkitmol2graph(mol_rdkit)
		fake_mol = FakeMolecule(mol_rdkit, mol_graph)
		self.assertTrue(len(get_unique_coords(fake_mol)) == 4)
		self.assertTrue(len(get_unique_coords(fake_mol, num_points=2)) == 3)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(0, 5))) == 3)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(1, 6))) == 4)
		for val, ind in zip(range(3), get_unique_coords(fake_mol, positions_to_check=range(0, 5))):
			self.assertEqual(val, ind[0])
		for val, ind in zip([1, 2, 4, 5], get_unique_coords(fake_mol, positions_to_check=range(1, 6))):
			self.assertEqual(val, ind[0])

		# non-ring
		mol_rdkit = Chem.MolFromSmiles('CC(N)C')
		mol_graph = rdkitmol2graph(mol_rdkit)
		fake_mol = FakeMolecule(mol_rdkit, mol_graph)
		self.assertTrue(len(get_unique_coords(fake_mol)) == 3)
		self.assertTrue(len(get_unique_coords(fake_mol, num_points=2)) == 0)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(0, 3))) == 3)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=range(1, 4))) == 3)
		self.assertTrue(len(get_unique_coords(fake_mol, positions_to_check=[0, 1, 4])) == 2)
		for val, ind in zip(range(3), get_unique_coords(fake_mol, positions_to_check=range(3))):
			self.assertEqual(val, ind[0])
		for val, ind in zip(range(1, 4), get_unique_coords(fake_mol, positions_to_check=range(1, 4))):
			self.assertEqual(val, ind[0])
		for val, ind in zip([0, 1], get_unique_coords(fake_mol, positions_to_check=[0, 1, 4])):
			self.assertEqual(val, ind[0])

	def test_3_unique_mols(self):
		mols_smiles = ['c1ccccc1', 'C1=CC=CC=C1', 'C1CCNCC1', 'C1NCCCC1', 'CC(N)C', 'NC(C)C']
		mols = []
		for s in mols_smiles:
			mol_rdkit = Chem.MolFromSmiles(s)
			mol_graph = rdkitmol2graph(mol_rdkit)
			mols.append(FakeModifier(mol=FakeMolecule(mol_rdkit, mol_graph)))
		self.assertTrue(len(get_unique_mols(mols)) == 3)
		mols_smiles = ['N1CCCCC1', 'C1NCCCC1', 'C1CNCCC1', 'C1CCNCC1', 'C1CCCNC1', 'C1CCCCN1']
		mols = []
		for s in mols_smiles:
			mol_rdkit = Chem.MolFromSmiles(s)
			mol_graph = rdkitmol2graph(mol_rdkit)
			mols.append(FakeModifier(mol=FakeMolecule(mol_rdkit, mol_graph)))
		self.assertTrue(len(get_unique_mols(mols)) == 1)


if __name__ == '__main__':
	unittest.main()
