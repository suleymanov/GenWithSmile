"""
Module for linkers generation between source and target molecules.
"""

class SourceTargetClient(object):
	def __init__(self, config):
		raise NotImplementedError()

	def process(self):
		raise NotImplementedError()

	@staticmethod
	def from_file(fn):
		raise NotImplementedError()
