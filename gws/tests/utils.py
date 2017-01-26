import os
import tempfile
import unittest
from shutil import rmtree
from json import dumps, loads
from random import sample

from rdkit import Chem

from gws.src.io import IOUtils


class ClientBaseTestCase(unittest.TestCase):
	def setUp(self):
		self.num_repeats = 2
		self.files_to_remove = []
		self.paths_to_remove = []
		self.atoms = IOUtils.read_smi('gws/tests/resources/atoms.smi')
		self.basic = IOUtils.read_smi('gws/tests/resources/basic.smi')
		self.groups = IOUtils.read_smi('gws/tests/resources/groups.smi')
		self.rings = IOUtils.read_smi('gws/tests/resources/rings.smi')

	def tearDown(self):
		# clean_files(self.files_to_remove)
		# clean_paths(self.paths_to_remove)
		pass


def clean_files(file_names):
	map(os.remove, filter(lambda fn: os.path.exists(fn), file_names))


def clean_paths(path_names):
	map(rmtree, filter(lambda pn: os.path.exists(pn), path_names))
