"""
Module for various combinations generation.
"""

class CombinationsClient(object):
	def __init__(self, config):
		raise NotImplementedError()

	def process(self):
		raise NotImplementedError()

	@staticmethod
	def from_file(fn):
		raise NotImplementedError()
