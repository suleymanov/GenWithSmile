class ValidationFactory(object):
	@staticmethod
	def validate_molecule(mol):
		assert isinstance(mol.smiles, str) or isinstance(mol.smiles, unicode), \
			'Molecule should be of type str.'
		
		assert isinstance(mol.atoms, list), 'Atoms should be provided as list: {}'.format(mol.atoms)
		assert isinstance(mol.attach_pos, list), \
			'Attach positions should be provided as list: {}'.format(mol.attach_pos)
		assert isinstance(mol.merge_pos, list), \
			'Merge positions should be provided as list: {}'.format(mol.merge_pos)
		assert all(map(lambda x: isinstance(x, str) or isinstance(x, unicode), mol.atoms)), \
			'Atom symbols should be of type str.'
		assert all(map(lambda x: isinstance(x, int), mol.attach_pos + mol.merge_pos)), \
			'Modification positions should be of type int.'

		assert mol.smiles != '', 'Molecule cannot be empty.'
		assert all(x != '' for x in mol.atoms), 'Atom symbols cannot be empty.'
		assert all(map(lambda x: x >= 0, mol.attach_pos + mol.merge_pos)), \
			'Modification positions cannot be negative.'
		assert len(mol.attach_pos) == len(set(mol.attach_pos)), \
			'Attach positions should be unique.'
		assert len(mol.merge_pos) == len(set(mol.merge_pos)), \
			'Merge positions should be unique.'

		from rdkit import Chem
		_tmp_mol = Chem.MolFromSmiles(mol.smiles)
		assert _tmp_mol is not None, 'Could not initialize molecule from SMILES: {}'.format(mol.smiles)
		assert all(map(lambda x: x < _tmp_mol.GetNumAtoms(), mol.attach_pos + mol.merge_pos)), \
			'Modification positions exceed number of atoms.'

	@staticmethod
	def validate_patterns(patterns):
		from rdkit import Chem
		for patt in patterns:
			assert Chem.MolFromSmarts(patt), \
				'Could not initialize pattern: {}'.format(patt)

	@staticmethod
	def validate_bound(bound):
		assert bound.low <= bound.high, 'Lower bound cannot be greater than higher bound.'

	@staticmethod
	def validate_iterations(iterations):
		assert isinstance(iterations, list), 'Iterations should be provided as list.'
		assert len(iterations) > 0, 'Should be at least one iteration.'
		for i, it in enumerate(iterations):
			ValidationFactory.validate_modification(it.attach)
			ValidationFactory.validate_modification(it.merge)
			assert isinstance(it.max, int), 'Number of modifications should be of type int.'
			assert it.max > 0, 'Number of modifications should be positive.'
			assert any([it.attach.one_point, it.attach.two_point, it.merge.one_point, it.merge.two_point]), \
				'No modification options were set in iteration {}.'.format(i)
			assert len(it.attach.addons) + len(it.merge.addons) > 0, \
				'No addons were set in iteration {}'.format(i)

	@staticmethod
	def validate_modification(mod):
		map(ValidationFactory.validate_molecule, mod.addons)
		assert isinstance(mod.one_point, bool), 'One-point option should be of type bool.'
		assert isinstance(mod.two_point, bool), 'Two-point option should be of type bool.'

	@staticmethod
	def validate_output(output):
		assert isinstance(output.path, str) or isinstance(output.path, unicode), \
			'Output path name should be of type str.'
		assert isinstance(output.alias, str) or isinstance(output.alias, unicode), \
			'Output file alias should be of type str.'
		assert isinstance(output.max, int), 'Maximum number of entries per file should be of type int.'

		assert output.path != '', 'Output path name cannot be empty.'
		assert output.alias != '', 'Output file alias cannot be empty.'
		assert output.max > 0, 'Maximum number of entries should be positive.'

	@staticmethod
	def validate_threading(num_threads):
		assert isinstance(num_threads, int), 'Number of threads should be of type int.'
		assert num_threads > 0, 'Should be at least one thread specified.'

		from multiprocessing import cpu_count

		count = cpu_count()
		assert num_threads <= count, 'Number of threads cannot be greater than {}.'.format(count)
