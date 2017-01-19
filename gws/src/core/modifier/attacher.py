from rdkit import Chem

from _base import BaseModifier, bonds_def


class Attacher(BaseModifier):
	def _combine_one_point(self, other_modifier, coord, other_coord):
		results = []
		host_pos, att_pos = coord[0], other_coord[0]
		host_atom_impl_val = self.mol_rdkit.GetAtomWithIdx(host_pos).GetImplicitValence()
		att_atom_impl_val = other_modifier.mol_rdkit.GetAtomWithIdx(att_pos).GetImplicitValence()
		for b_order in bonds_def:
			if b_order <= host_atom_impl_val and b_order <= att_atom_impl_val:
				em = Chem.EditableMol(self.mol_rdkit)
				inds_map = {}
				for at in other_modifier.mol_rdkit.GetAtoms():
					inds_map[at.GetIdx()] = em.AddAtom(at)
				for b in other_modifier.mol_rdkit.GetBonds():
					em.AddBond(
						inds_map[b.GetBeginAtomIdx()], 
						inds_map[b.GetEndAtomIdx()], 
						b.GetBondType())
				em.AddBond(host_pos, inds_map[att_pos], bonds_def[b_order])
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

	def _combine_two_points(self, other_modifier, coord, other_coord):
		results = []
		host_impl_val = [
			self.mol_rdkit.GetAtomWithIdx(coord[0]).GetImplicitValence(),
			self.mol_rdkit.GetAtomWithIdx(coord[1]).GetImplicitValence()]
		att_impl_val = [
			other_modifier.mol_rdkit.GetAtomWithIdx(other_coord[0]).GetImplicitValence(),
			other_modifier.mol_rdkit.GetAtomWithIdx(other_coord[1]).GetImplicitValence()]
		for b_order in bonds_def:
			if b_order <= host_impl_val[0] and b_order <= att_impl_val[0]:
				for b_order2 in bonds_def:
					if b_order2 <= host_impl_val[1] and b_order2 <= att_impl_val[1]:
						em = Chem.EditableMol(self.mol_rdkit)
						inds_map = {}
						for at in other_modifier.mol_rdkit.GetAtoms():
							inds_map[at.GetIdx()] = em.AddAtom(at)
						for b in other_modifier.mol_rdkit.GetBonds():
							em.AddBond(
								inds_map[b.GetBeginAtomIdx()],
								inds_map[b.GetEndAtomIdx()],
								b.GetBondType())
						em.AddBond(coord[0], inds_map[other_coord[0]], bonds_def[b_order])
						em.AddBond(coord[1], inds_map[other_coord[1]], bonds_def[b_order2])
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
