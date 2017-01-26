import networkx as nx
from rdkit import Chem

from algs import rdkitmol2graph
from algs import reduce_graph
from algs import edge_type_to_label
from modifier import Attacher, Merger


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
		self.mol_smiles = _mol_smiles
		self.mol_rdkit = _mol_rdkit
		self.mol_graph = None
		self.attacher = Attacher(_mol_rdkit, attach_pos, atoms)
		self.merger = Merger(_mol_rdkit, merge_pos, atoms)

	def attach(self, addon_handler, one_point=True, two_point=False):
		"""
		Perform all unique attaches of some addon.
		:param addon_handler: MoleculeHandler
		:param alias: str
		:param one_point: bool
		:param two_point: bool
		:return: list of MoleculeHandler
		"""
		attach_results = []
		for res in self.attacher.combine(addon_handler.attacher, one_point=one_point, two_point=False):
			next_attach_pos = res.next_pos
			next_merge_pos = res.next_pos
			handler = MoleculeHandler._create_handler(
				res.result_mol_rdkit, res.pos, self.merger.pos, next_attach_pos, next_merge_pos)
			coord = self.attacher.get_coords(one_point=one_point)[res.ind]
			handler._update_mol_graph(self._create_mol_graph_from_mod(
				addon_handler.mol_graph, coord, res.result_mol_rdkit, res.inds_map))
			attach_results.append(handler)
		for res in self.attacher.combine(addon_handler.attacher, one_point=False, two_point=two_point):
			next_attach_pos = res.next_pos
			next_merge_pos = res.next_pos
			handler = MoleculeHandler._create_handler(
				res.result_mol_rdkit, res.pos, self.merger.pos, next_attach_pos, next_merge_pos)
			coord = self.attacher.get_coords(two_point=two_point)[res.ind]
			handler._update_mol_graph(self._create_mol_graph_from_mod(
				addon_handler.mol_graph, coord, res.result_mol_rdkit, res.inds_map))
			attach_results.append(handler)
		return attach_results

	def merge(self, addon_handler, one_point=True, two_point=False):
		"""
		Perform all unique merges of some addon.
		:param addon_handler: MoleculeHandler
		:param alias: str
		:param one_point: bool
		:param two_point: bool
		:return: list of MoleculeHandler
		"""
		merge_results = []
		for res in self.merger.combine(addon_handler.merger, one_point=one_point, two_point=False):
			next_attach_pos = res.next_pos
			next_merge_pos = res.next_pos
			handler = MoleculeHandler._create_handler(
				res.result_mol_rdkit, self.attacher.pos, res.pos, next_attach_pos, next_merge_pos)
			coord = self.merger.get_coords(one_point=one_point)[res.ind]
			handler._update_mol_graph(self._create_mol_graph_from_mod(
				addon_handler.mol_graph, coord, res.result_mol_rdkit, res.inds_map))
			merge_results.append(handler)
		for res in self.merger.combine(addon_handler.merger, one_point=False, two_point=two_point):
			next_attach_pos = res.next_pos
			next_merge_pos = res.next_pos
			handler = MoleculeHandler._create_handler(
				res.result_mol_rdkit, self.attacher.pos, res.pos, next_attach_pos, next_merge_pos)
			coord = self.merger.get_coords(two_point=two_point)[res.ind]
			handler._update_mol_graph(self._create_mol_graph_from_mod(
				addon_handler.mol_graph, coord, res.result_mol_rdkit, res.inds_map))
			merge_results.append(handler)
		return merge_results

	def update_modifiers_positions(self):
		self.attacher.update_pos()
		self.merger.update_pos()

	def _update_mol_graph(self, mol_graph):
		self.mol_graph = mol_graph
		self.attacher.set_mol_graph(mol_graph)
		self.merger.set_mol_graph(mol_graph)

	def create_mol_graph(self, alias=None):
		"""
		Create reduced molecular graph from rdkit molecule.
		:param alias: str
		:return: None
		"""
		self.mol_graph = reduce_graph(
			rdkitmol2graph(self.mol_rdkit), list(set(self.attacher.pos + self.merger.pos)), alias)
		# self.mol_graph = rdkitmol2graph(self.mol_rdkit)
		self._update_mol_graph(self.mol_graph)

	def _create_mol_graph_from_mod(self, other_mol_graph, coord, result_mol, inds_map):
		"""
		Create reduced molecular graph from modification result.
		:param other_mol_graph: nx.Graph
		:param coord: list of int
		:param other_coord: list of int
		:param result_mol: Chem.rdchem.Mol
		:param inds_map: dict
		:return: nx.Graph
		"""
		assert all(map(lambda x: x in self.mol_graph.nodes(), coord))

		new_graph = nx.Graph(self.mol_graph.copy())

		update_node = lambda node: inds_map[node] if node in inds_map else node

		for node in filter(lambda x: x in inds_map, other_mol_graph.nodes()):
			new_graph.add_node(inds_map[node], label=other_mol_graph.node[node]['label'])
		for node in filter(lambda x: x not in inds_map, other_mol_graph.nodes()):
			new_graph.add_node(
				'{}_other'.format(node) if node in new_graph.nodes() else node, 
				label=other_mol_graph.node[node]['label'])
		for edge in other_mol_graph.edges():
			edge_data = other_mol_graph.get_edge_data(*edge)
			new_graph.add_edge(
				*map(update_node, edge), weight=edge_data['weight'], label=edge_data['label'])
		for c in coord:
			for node in filter(lambda x: x in inds_map, other_mol_graph.nodes()):
				bond = result_mol.GetBondBetweenAtoms(c, inds_map[node])
				if bond:
					edge_type = int(bond.GetBondType())
					label = edge_type_to_label.get(edge_type, '')
					new_graph.add_edge(c, inds_map[node], weight=edge_type, label=label)
		return new_graph

	@staticmethod
	def _create_handler(mol_rdkit, attach_pos, merge_pos, next_attach_pos, next_merge_pos):
		handler = MoleculeHandler(mol_rdkit=mol_rdkit, attach_pos=attach_pos, merge_pos=merge_pos)
		handler.attacher.add_next_pos(next_attach_pos)
		handler.merger.add_next_pos(next_merge_pos)
		return handler
