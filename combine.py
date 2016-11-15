from itertools import izip

from rdkit import Chem

# from utils import rdkitmol2graph


bonds_def = {
	1.0: Chem.rdchem.BondType.SINGLE,
	2.0: Chem.rdchem.BondType.DOUBLE,
	3.0: Chem.rdchem.BondType.TRIPLE,
	4.0: Chem.rdchem.BondType.QUADRUPLE,
	5.0: Chem.rdchem.BondType.QUINTUPLE
}


def combine_mols(mol_smiles_1, mol_smiles_2, pos_1, pos_2, src_atom, trg_atom):
	"""
	Combine two molecules provided modification positions.
	:param mol_smiles_1: str
	:param mol_smiles_2: str
	:param pos_1: int
	:param pos_2: int
	:param src_atom: str
	:param trg_atom: str
	:return: list of str
	"""
	combination_results = []
	mol_1 = Chem.MolFromSmiles(mol_smiles_1)
	mol_2 = Chem.MolFromSmiles(mol_smiles_2)
	if (mol_1.GetAtomWithIdx(pos_1).GetSymbol() == src_atom 
		and mol_2.GetAtomWithIdx(pos_2).GetSymbol() == trg_atom):
		Chem.Kekulize(mol_1)
		Chem.Kekulize(mol_2)
		combination_results.extend(_perform_attach(mol_1, mol_2, pos_1, pos_2))
	return combination_results


def merge_bond(mol_smiles_1, mol_smiles_2, bond_1, bond_2):
	"""
	Combine two molecules by merging specified bond into one in all possible ways.
	:param mol_smiles_1: str
	:param mol_smiles_2: str
	:param bond_1: tuple of int
	:param bond_2: tuple of int
	:return: list of strs
	"""
	combination_results = []
	mol_1 = Chem.MolFromSmiles(mol_smiles_1)
	mol_2 = Chem.MolFromSmiles(mol_smiles_2)
	Chem.Kekulize(mol_1)
	Chem.Kekulize(mol_2)
	combination_results.extend(_merge_bond(mol_1, mol_2, bond_1, bond_2))
	combination_results.extend(_merge_bond(mol_1, mol_2, bond_1, bond_2, turn=True))
	return combination_results


def _merge_bond(host_mol, merged_mol, h_bonds, m_bonds, turn=False):
	def get_bond(mol_bonds, bond_beg, bond_end):
		bonds = filter(
			lambda b: b.GetBeginAtomIdx() == bond_beg and b.GetEndAtomIdx() == bond_end, mol_bonds)
		assert len(bonds) == 1
		return bonds[0]

	def is_acceptable_bond_order(
		_b_order, _h_b_at_impl_val, _h_e_at_impl_val, _rm_b_at_exp_val, _rm_e_at_exp_val, _turn):
		if _turn:
			return (
				_h_b_at_impl_val >= _rm_e_at_exp_val + _b_order
				and _h_e_at_impl_val >= _rm_b_at_exp_val + _b_order
			)
		else:
			return (
				_h_b_at_impl_val >= _rm_b_at_exp_val + _b_order
				and _h_e_at_impl_val >= _rm_e_at_exp_val + _b_order
			)

	results = []
	host_dval = get_bond(host_mol.GetBonds(), *h_bonds).GetBondTypeAsDouble()
	rem_dval = get_bond(merged_mol.GetBonds(), *m_bonds).GetBondTypeAsDouble()

	h_b_at_impl_val = host_mol.GetAtomWithIdx(h_bonds[0]).GetImplicitValence() + host_dval
	h_e_at_impl_val = host_mol.GetAtomWithIdx(h_bonds[1]).GetImplicitValence() + host_dval
	rm_b_at_exp_val = merged_mol.GetAtomWithIdx(m_bonds[0]).GetExplicitValence() - rem_dval
	rm_e_at_exp_val = merged_mol.GetAtomWithIdx(m_bonds[1]).GetExplicitValence() - rem_dval
	
	for b_order in bonds_def:
		if is_acceptable_bond_order(
			b_order, h_b_at_impl_val, h_e_at_impl_val, rm_b_at_exp_val, rm_e_at_exp_val, turn):
			em = Chem.EditableMol(host_mol)
			em.RemoveBond(*h_bonds)
			new_atoms_inds_map = (
				{m_bonds[0]: h_bonds[1], m_bonds[1]: h_bonds[0]} if turn else 
				{m_bonds[0]: h_bonds[0], m_bonds[1]: h_bonds[1]}
			)  # atom: old position -> new position
			for at in filter(lambda at: at.GetIdx() not in m_bonds, merged_mol.GetAtoms()):
				new_atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
			for b in filter(
				lambda b: b.GetBeginAtomIdx() != m_bonds[0] and b.GetEndAtomIdx() != m_bonds[1], 
				merged_mol.GetBonds()):
				em.AddBond(
					new_atoms_inds_map[b.GetBeginAtomIdx()],
					new_atoms_inds_map[b.GetEndAtomIdx()],
					b.GetBondType()
				)
			em.AddBond(h_bonds[0], h_bonds[1], bonds_def[b_order])
			result_mol = em.GetMol()
			Chem.SanitizeMol(result_mol)
			results.append(Chem.MolToSmiles(result_mol))
	return results


def _perform_insertion(host_mol, ins_mol, host_pos, ins_pos):
	results = []
	host_atom_exp_val = host_mol.GetAtomWithIdx(host_pos).GetExplicitValence()
	ins_atom_imp_val = ins_mol.GetAtomWithIdx(ins_pos).GetImplicitValence()
	if host_atom_exp_val <= ins_atom_imp_val:
		em = Chem.EditableMol(host_mol)
		new_atoms_inds_map = {ins_pos: host_pos}  # atom: old position -> new position
		em.ReplaceAtom(host_pos, ins_mol.GetAtomWithIdx(ins_pos))
		for at in filter(lambda x: x.GetIdx() != ins_pos, ins_mol.GetAtoms()):
			new_atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
		for b in ins_mol.GetBonds():
			em.AddBond(
				new_atoms_inds_map[b.GetBeginAtomIdx()],
				new_atoms_inds_map[b.GetEndAtomIdx()],
				b.GetBondType()
			)
		result_mol = em.GetMol()
		if host_mol.GetAtomWithIdx(host_pos).GetIsAromatic():
			result_mol.GetAtomWithIdx(host_pos).SetIsAromatic(True)
		Chem.SanitizeMol(result_mol)
		results.append(Chem.MolToSmiles(result_mol))
	return results


def _perform_attach(host_mol, att_mol, host_pos, att_pos):
	results = []
	host_atom_impl_val = host_mol.GetAtomWithIdx(host_pos).GetImplicitValence()
	att_atom_impl_val = att_mol.GetAtomWithIdx(att_pos).GetImplicitValence()
	for b_order in bonds_def:
		if b_order <= host_atom_impl_val and b_order <= att_atom_impl_val:
			em = Chem.EditableMol(host_mol)
			new_atoms_inds_map = {}  # atom: old position -> new position
			for at in att_mol.GetAtoms():
				new_atoms_inds_map[at.GetIdx()] = em.AddAtom(at)
			for b in att_mol.GetBonds():
				em.AddBond(
					new_atoms_inds_map[b.GetBeginAtomIdx()],
					new_atoms_inds_map[b.GetEndAtomIdx()],
					b.GetBondType()
				)
			em.AddBond(host_pos, new_atoms_inds_map[att_pos], bonds_def[b_order])
			result_mol = em.GetMol()
			Chem.SanitizeMol(result_mol)
			results.append(Chem.MolToSmiles(result_mol))
	return results
