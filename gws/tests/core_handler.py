import unittest

from rdkit import Chem

from gws.src.core import MoleculeHandler
from gws.src.core.algs import rdkitmol2graph


class CoreHandlerTests(unittest.TestCase):
	def test_ctor_1(self):
		with self.assertRaises(StandardError):
			MoleculeHandler()

		mol_smiles = 'c1ccccc1'
		MoleculeHandler(mol_smiles=mol_smiles)
		MoleculeHandler(mol_rdkit=Chem.MolFromSmiles(mol_smiles))

	def test_attach_1(self):
		"""
		Single attach in one point (cyclohexane with another cyclohexane).
		"""
		handler = MoleculeHandler('C1CCCCC1')
		results = handler.attach(handler)
		self.assertTrue(len(results) == 2)
		for res in results:
			self.assertTrue(len(res.attacher.pos) == 5)
			self.assertTrue(len(res.attacher.next_pos) == 6)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 6)
			res.update_modifiers_positions()
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 0)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 0)

	def test_attach_2(self):
		"""
		Single attach in two points (cyclohexane with another cyclohexane).
		"""
		handler = MoleculeHandler('C1CCCCC1')
		results = handler.attach(handler, one_point=False, two_point=True)
		self.assertTrue(len(results) == 4)
		for res in results:
			self.assertTrue(len(res.attacher.pos) == 4)
			self.assertTrue(len(res.attacher.next_pos) == 6)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 6)
			res.update_modifiers_positions()
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 0)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 0)

	def test_merge_1(self):
		"""
		Single merge in one point (cyclohexane to another cyclohexane).
		"""
		handler = MoleculeHandler('C1CCCCC1')
		results = handler.merge(handler)
		self.assertTrue(len(results) == 1)
		for res in results:
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 6)
			self.assertTrue(len(res.merger.pos) == 5)
			self.assertTrue(len(res.merger.next_pos) == 6)
			res.update_modifiers_positions()
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 0)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 0)

	def test_merge_2(self):
		"""
		Single merge in two points (cyclohexane to another cyclohexane).
		"""
		handler = MoleculeHandler('C1CCCCC1')
		results = handler.merge(handler, one_point=False, two_point=True)
		self.assertTrue(len(results) == 2)
		for res in results:
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 6)
			self.assertTrue(len(res.merger.pos) == 4)
			self.assertTrue(len(res.merger.next_pos) == 6)
			res.update_modifiers_positions()
			self.assertTrue(len(res.attacher.pos) == 6)
			self.assertTrue(len(res.attacher.next_pos) == 0)
			self.assertTrue(len(res.merger.pos) == 6)
			self.assertTrue(len(res.merger.next_pos) == 0)


if __name__ == '__main__':
	unittest.main()
