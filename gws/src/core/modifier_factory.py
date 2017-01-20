"""
Module for modifying molecule.
"""

from rdkit import Chem


bonds_def = {
	1.0: Chem.rdchem.BondType.SINGLE,
	2.0: Chem.rdchem.BondType.DOUBLE,
	3.0: Chem.rdchem.BondType.TRIPLE,
	4.0: Chem.rdchem.BondType.QUADRUPLE,
	5.0: Chem.rdchem.BondType.QUINTUPLE
}


class ModifierFactory(object):
	@staticmethod
	def merge(host, addon, host_coord, addon_coord):
		"""
		Generate all possible merges of 'addon.mol' into 'host.mol' with coordinates provided.
		Results positions are updated appropriately.
		:param host: modifier.Modifier
		:param addon: modifier.Modifier
		:param host_coord: list of int (length == 1 or 2)
		:param addon_coord: list of int (length == 1 or 2)
		:return: list of modifier.Modifier
		"""
		assert len(host_coord) == len(addon_coord)
		if len(host_coord) == 1:
			return ModifierFactory._merge_atom(host, addon, host_coord, addon_coord)
		elif len(host_coord) == 2:
			return ModifierFactory._merge_bond(host, addon, host_coord, addon_coord)
		raise ValueError('Three and more point modifications await their implementation.')

	@staticmethod
	def attach(host, addon, host_coord, addon_coord):
		"""
		Generate all possible attaches of 'addon.mol' into 'host.mol' with coordinates provided.
		Results positions are updated appropriately.
		:param host: modifier.Modifier
		:param addon: modifier.Modifier
		:param host_coord: list of int (length == 1 or 2)
		:param addon_coord: list of int (length == 1 or 2)
		:return: list of modifier.Modifier
		"""
		assert len(host_coord) == len(addon_coord)
		if len(host_coord) == 1:
			return ModifierFactory._attach(host, addon, host_coord, addon_coord)
		if len(host_coord) == 2:
			return ModifierFactory._two_point_attach(host, addon, host_coord, addon_coord)
		raise ValueError('Three and more point modifications await their implementation.')

	@staticmethod
	def _attach(host_modifier, addon_modifier, host_coord, addon_coord):
		host_pos, att_pos = host_coord[0], addon_coord[0]
		results = []
		host_atom_impl_val = host_modifier.mol.rdkit.GetAtomWithIdx(host_pos).GetImplicitValence()
		att_atom_impl_val = addon_modifier.mol.rdkit.GetAtomWithIdx(att_pos).GetImplicitValence()
		for b_order in bonds_def:
			if b_order <= host_atom_impl_val and b_order <= att_atom_impl_val:
				em = Chem.EditableMol(host_modifier.mol.rdkit)
				atoms_inds_map = {}  # old index : new index
				for at in addon_modifier.mol.rdkit.GetAtoms():
					atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
				for b in addon_modifier.mol.rdkit.GetBonds():
					em.AddBond(
						atoms_inds_map[b.GetBeginAtomIdx()],
						atoms_inds_map[b.GetEndAtomIdx()],
						b.GetBondType())
				em.AddBond(host_pos, atoms_inds_map[att_pos], bonds_def[b_order])
				result_mol = em.GetMol()
				try:
					Chem.SanitizeMol(result_mol)
				except ValueError as e:
					print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
					print(e.message)
				except StandardError as e:
					print(e.message)
				else:
					pop_pos = atoms_inds_map.pop(att_pos)
					results.append(host_modifier.create_new_modifier(
						result_mol, addon_modifier, atoms_inds_map, rm_attach_pos=host_coord,
						pop_pos=pop_pos))
		return results

	@staticmethod
	def _two_point_attach(host_modifier, addon_modifier, host_coord, addon_coord):
		results = []
		host_impl_val = [
			host_modifier.mol.rdkit.GetAtomWithIdx(host_coord[0]).GetImplicitValence(),
			host_modifier.mol.rdkit.GetAtomWithIdx(host_coord[1]).GetImplicitValence()]
		att_impl_val = [
			addon_modifier.mol.rdkit.GetAtomWithIdx(addon_coord[0]).GetImplicitValence(),
			addon_modifier.mol.rdkit.GetAtomWithIdx(addon_coord[1]).GetImplicitValence()]
		for b_order in bonds_def:
			if b_order <= host_impl_val[0] and b_order <= att_impl_val[0]:
				for b_order2 in bonds_def:
					if b_order2 <= host_impl_val[1] and b_order2 <= att_impl_val[1]:
						em = Chem.EditableMol(host_modifier.mol.rdkit)
						atoms_inds_map = {}  # old index : new index
						for at in addon_modifier.mol.rdkit.GetAtoms():
							atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
						for b in addon_modifier.mol.rdkit.GetBonds():
							em.AddBond(
								atoms_inds_map[b.GetBeginAtomIdx()],
								atoms_inds_map[b.GetEndAtomIdx()],
								b.GetBondType())
						em.AddBond(host_coord[0], atoms_inds_map[addon_coord[0]], bonds_def[b_order])
						em.AddBond(host_coord[1], atoms_inds_map[addon_coord[1]], bonds_def[b_order2])
						result_mol = em.GetMol()
						try:
							Chem.SanitizeMol(result_mol)
						except ValueError as e:
							print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
							print(e.message)
						except StandardError as e:
							print(e.message)
						else:
							results.append(host_modifier.create_new_modifier(
								result_mol, addon_modifier, atoms_inds_map, rm_attach_pos=host_coord))
		return results

	@staticmethod
	def _merge_atom(host_modifier, addon_modifier, host_coord, addon_coord):
		host_pos, ins_pos = host_coord[0], addon_coord[0]
		results = []
		host_atom_exp_val = host_modifier.mol.rdkit.GetAtomWithIdx(host_pos).GetExplicitValence()
		ins_atom_imp_val = addon_modifier.mol.rdkit.GetAtomWithIdx(ins_pos).GetImplicitValence()
		if host_atom_exp_val <= ins_atom_imp_val:
			em = Chem.EditableMol(host_modifier.mol.rdkit)
			atoms_inds_map = {ins_pos: host_pos}  # old index : new index
			em.ReplaceAtom(host_pos, addon_modifier.mol.rdkit.GetAtomWithIdx(ins_pos))
			for at in filter(lambda x: x.GetIdx() != ins_pos, addon_modifier.mol.rdkit.GetAtoms()):
				atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
			for b in addon_modifier.mol.rdkit.GetBonds():
				em.AddBond(
					atoms_inds_map[b.GetBeginAtomIdx()],
					atoms_inds_map[b.GetEndAtomIdx()],
					b.GetBondType())
			result_mol = em.GetMol()
			if host_modifier.mol.rdkit.GetAtomWithIdx(host_pos).GetIsAromatic():
				result_mol.GetAtomWithIdx(host_pos).SetIsAromatic(True)
			try:
				Chem.SanitizeMol(result_mol)
			except ValueError as e:
				print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
				print(e.message)
			except StandardError as e:
				print(e.message)
			else:
				pop_pos = atoms_inds_map.pop(ins_pos)
				results.append(host_modifier.create_new_modifier(
					result_mol, addon_modifier, atoms_inds_map, rm_merge_pos=host_coord,
					pop_pos=pop_pos))
		return results

	@staticmethod
	def _merge_bond(host_modifier, addon_modifier, host_coord, addon_coord):
		def get_bond(mol_bonds, bond_beg, bond_end):
			bonds = filter(
				lambda b: b.GetBeginAtomIdx() == bond_beg and b.GetEndAtomIdx() == bond_end, 
				mol_bonds)
			assert len(bonds) == 1
			return bonds[0]

		def is_acceptable_bond_order(
			_b_order, _h_b_at_impl_val, _h_e_at_impl_val, _rm_b_at_exp_val, _rm_e_at_exp_val):
			return (
				_h_b_at_impl_val >= _rm_b_at_exp_val + _b_order and
				_h_e_at_impl_val >= _rm_e_at_exp_val + _b_order)

		results = []
		host_dval = get_bond(host_modifier.mol.rdkit.GetBonds(), *host_coord).GetBondTypeAsDouble()
		rem_dval = get_bond(addon_modifier.mol.rdkit.GetBonds(), *addon_coord).GetBondTypeAsDouble()
		h_b_at_impl_val = host_modifier.mol.rdkit.GetAtomWithIdx(host_coord[0]).GetImplicitValence() + host_dval
		h_e_at_impl_val = host_modifier.mol.rdkit.GetAtomWithIdx(host_coord[1]).GetImplicitValence() + host_dval
		rm_b_at_exp_val = addon_modifier.mol.rdkit.GetAtomWithIdx(addon_coord[0]).GetExplicitValence() - rem_dval
		rm_e_at_exp_val = addon_modifier.mol.rdkit.GetAtomWithIdx(addon_coord[1]).GetExplicitValence() - rem_dval
		for b_order in bonds_def:
			if is_acceptable_bond_order(
				b_order, h_b_at_impl_val, h_e_at_impl_val, rm_b_at_exp_val, rm_e_at_exp_val):
				em = Chem.EditableMol(host_modifier.mol.rdkit)
				em.RemoveBond(*host_coord)
				atoms_inds_map = {addon_coord[i]: host_coord[i] for i in [0, 1]}  # old : new
				for at in filter(lambda at: at.GetIdx() not in addon_coord, addon_modifier.mol.rdkit.GetAtoms()):
					atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
				for b in filter(
					lambda b: b.GetBeginAtomIdx() != addon_coord[0] and b.GetEndAtomIdx() != addon_coord[1],
					addon_modifier.mol.rdkit.GetBonds()):
					em.AddBond(
						atoms_inds_map[b.GetBeginAtomIdx()],
						atoms_inds_map[b.GetEndAtomIdx()],
						b.GetBondType())
				em.AddBond(host_coord[0], host_coord[1], bonds_def[b_order])
				result_mol = em.GetMol()
				try:
					Chem.SanitizeMol(result_mol)
				except ValueError as e:
					print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
					print(e.message)
				except StandardError as e:
					print(e.message)
				else:
					results.append(host_modifier.create_new_modifier(
						result_mol, addon_modifier, atoms_inds_map, rm_merge_pos=host_coord))
		return results
