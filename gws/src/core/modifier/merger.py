from rdkit import Chem

from _base import BaseModifier, bonds_def


class Merger(BaseModifier):
	def _combine_one_point(self, other_modifier, coord, other_coord):
		results = []
		host_pos, ins_pos = coord[0], other_coord[0]
		host_atom_exp_val = self.mol_rdkit.GetAtomWithIdx(host_pos).GetExplicitValence()
		ins_atom_imp_val = other_modifier.mol_rdkit.GetAtomWithIdx(ins_pos).GetImplicitValence()
		if host_atom_exp_val <= ins_atom_imp_val:
			em = Chem.EditableMol(self.mol_rdkit)
			inds_map = {ins_pos: host_pos}
			em.ReplaceAtom(host_pos, other_modifier.mol_rdkit.GetAtomWithIdx(ins_pos))
			for at in filter(lambda x: x.GetIdx() != ins_pos, other_modifier.mol_rdkit.GetAtoms()):
				inds_map[at.GetIdx()] = em.AddAtom(at)
			for b in other_modifier.mol_rdkit.GetBonds():
				em.AddBond(
					inds_map[b.GetBeginAtomIdx()],
					inds_map[b.GetEndAtomIdx()],
					b.GetBondType())
			result_mol = em.GetMol()
			if self.mol_rdkit.GetAtomWithIdx(host_pos).GetIsAromatic():
				result_mol.GetAtomWithIdx(host_pos).SetIsAromatic(True)
			try:
				Chem.SanitizeMol(result_mol)
			except ValueError as e:
				print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
				print(e.message)
			except StandardError as e:
				print(e.message)
			else:
				results.append(self._pack_result(other_modifier, result_mol, inds_map, coord))
		return results

	def _combine_two_points(self, other_modifier, coord, other_coord):
		get_bond = lambda bonds, beg, end: filter(
			lambda x: x.GetBeginAtomIdx() == beg and x.GetEndAtomIdx() == end, bonds)[0]
		is_acceptable_bond_order = lambda b_order, h_b_impl, h_e_impl, rm_b_exp, rm_e_exp: (
			h_b_impl >= rm_b_exp + b_order and h_e_impl >= rm_e_exp + b_order)

		results = []
		dval = get_bond(self.mol_rdkit.GetBonds(), *coord).GetBondTypeAsDouble()
		rm_dval = get_bond(other_modifier.mol_rdkit.GetBonds(), *other_coord).GetBondTypeAsDouble()
		h_b_impl = self.mol_rdkit.GetAtomWithIdx(coord[0]).GetImplicitValence() + dval
		h_e_impl = self.mol_rdkit.GetAtomWithIdx(coord[1]).GetImplicitValence() + dval
		rm_b_exp = other_modifier.mol_rdkit.GetAtomWithIdx(other_coord[0]).GetImplicitValence() - rm_dval
		rm_e_exp = other_modifier.mol_rdkit.GetAtomWithIdx(other_coord[1]).GetImplicitValence() - rm_dval
		for b_order in bonds_def:
			if is_acceptable_bond_order(b_order, h_b_impl, h_e_impl, rm_b_exp, rm_e_exp):
				em = Chem.EditableMol(self.mol_rdkit)
				em.RemoveBond(*coord)
				inds_map = {other_coord[i]: coord[i] for i in [0, 1]}
				for at in filter(lambda at: at.GetIdx() not in other_coord, other_modifier.mol_rdkit.GetAtoms()):
					inds_map[at.GetIdx()] = em.AddAtom(at)
				for b in filter(
					lambda b: b.GetBeginAtomIdx() != other_coord[0] and b.GetEndAtomIdx() != other_coord[1],
					other_modifier.mol_rdkit.GetBonds()):
					em.AddBond(
						inds_map[b.GetBeginAtomIdx()],
						inds_map[b.GetEndAtomIdx()],
						b.GetBondType())
				em.AddBond(coord[0], coord[1], bonds_def[b_order])
				result_mol = em.GetMol()
				try:
					Chem.SanitizeMol(result_mol)
				except ValueError as e:
					print('Could not sanitize molecule: {}'.format(Chem.MolToSmiles(result_mol)))
					print(e.message)
				except StandardError as e:
					print(e.message)
				else:
					results.append(self._pack_result(other_modifier, result_mol, inds_map, coord))
		return results
