import os
from json import dumps, loads


class IOUtils(object):
	@staticmethod
	def read_file(fn):
		assert os.path.exists(fn)

		with open(fn) as f:
			return f.read()

	@staticmethod
	def write_file(fn, data):
		with open(fn, 'w') as f:
			f.write(data)

	@staticmethod
	def read_json(fn):
		return loads(IOUtils.read_file(fn))

	@staticmethod
	def read_smi(fn):
		return filter(lambda x: x != '', IOUtils.read_file(fn).splitlines())

	@staticmethod
	def write_json(fn, data):
		return IOUtils.write_file(fn, dumps(data, indent=4))

	@staticmethod
	def write_smi(fn, data):
		return IOUtils.write_file(fn, '\n'.join(data) + '\n')

	@staticmethod
	def read_df(fn, sep='|', col_ind=1):
		import pandas as pd

		df = pd.read_csv(fn, sep=sep)
		return list(df.iloc[:,col_ind])
