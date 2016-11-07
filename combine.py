from rdkit import Chem

from utils import rdkitmol2graph


bonds_def = {
	1: Chem.rdchem.BondType.SINGLE,
	2: Chem.rdchem.BondType.DOUBLE,
	3: Chem.rdchem.BondType.TRIPLE,
	4: Chem.rdchem.BondType.QUADRUPLE,
	5: Chem.rdchem.BondType.QUINTUPLE
}


def combine_mols(mol_smiles_1, mol_smiles_2, pos_1, pos_2):
	"""
	Combine two molecules provided modification positions.
	:param mol_smiles_1: str
	:param mol_smiles_2: str
	:param pos_1: int
	:param pos_2: int
	:return: list of str
	"""
	combination_results = []
	mol_1 = Chem.MolFromSmiles(mol_smiles_1)
	mol_2 = Chem.MolFromSmiles(mol_smiles_2)
	combination_results.extend(_perform_insertion(mol_1, mol_2, pos_1, pos_2))
	combination_results.extend(_perform_insertion(mol_2, mol_1, pos_2, pos_1))
	combination_results.extend(_perform_attach(mol_1, mol_2, pos_1, pos_2))
	return combination_results


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
