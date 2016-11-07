"""
Module for molecule representation
"""

from rdkit import Chem

from algs import rdkitmol2graph


class Molecule(object):
	def __init__(self, mol_smiles=None, mol_rdkit=None, mol_graph=None, positions=None):
		"""
		Create object. Assumes all inputs are fully compatible with each other (if not None).
		:param mol_smiles: str
		:param mol_rdkit: Chem.Mol
		:param mol_graph: networkx.Graph
		:param positions: list of int
		:return: None
		"""
		self.mol_smiles = mol_smiles
		self.mol_rdkit = mol_rdkit
		self.mol_graph = mol_graph
		self.positions = positions
		self._next_positions = []

	def add_next_positions(self, positions):
		self._next_positions = list(set(self._next_positions + positions))

	def update_positions(self):
		self.positions = self._next_positions[::]
		self._next_positions = []

	def get_smiles(self):
		return self.mol_smiles if self.mol_smiles else Chem.MolToSmiles(self.mol_rdkit)

	@staticmethod
	def from_smiles(mol_smiles, allowed_symbols=None):
		"""
		Initialize from molecule in SMILES format.
		:param mol_smiles: str
		:param allowed_symbols: list of str
		:return: Molecule
		"""
		mol_rdkit = Chem.MolFromSmiles(mol_smiles)
		mol_graph = rdkitmol2graph(mol_rdkit)
		positions = Molecule._define_positions(mol_rdkit, allowed_symbols)
		return Molecule(
			mol_smiles=mol_smiles, 
			mol_rdkit=mol_rdkit, 
			mol_graph=mol_graph, 
			positions=positions)

	@staticmethod
	def from_rdkit(mol_rdkit, positions=None):
		"""
		Initialize from molecule in rdkit format.
		:param mol_rdkit: Chem.Mol
		:return: Molecule
		"""
		# raise NotImplementedError()

		mol_graph = rdkitmol2graph(mol_rdkit)
		return Molecule(mol_rdkit=mol_rdkit, mol_graph=mol_graph, positions=positions)

	@staticmethod
	def from_graph(mol_graph):
		"""
		Initialize from molecule in graph format.
		:param mol_graph: networkx.Graph
		:return: Molecule
		"""
		raise NotImplementedError()

	@staticmethod
	def _define_positions(mol, allowed_symbols):
		return (
			map(
				lambda at: at.GetIdx(),
				filter(lambda at: at.GetSymbol() in allowed_symbols, mol.GetAtoms()))
			if allowed_symbols else list(xrange(mol.GetNumAtoms())))