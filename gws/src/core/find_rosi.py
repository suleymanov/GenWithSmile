import os
import sys

from rdkit import Chem

from algs import get_unique_mols, rdkitmol2graph, is_isomorph

s = 'O=C1NC(=O)SC1Cc3ccc(OCCN(c2ncccc2)C)cc3'
# res_pn = '../../../test_run/rosi_test/'
res_pn = 'tmp/'
assert os.path.exists(res_pn)


def read_smi(fn):
	with open(fn) as f:
		return f.read().splitlines()

results = map(
	lambda i: reduce(
		lambda res, item: res + read_smi(res_pn + item),
		filter(lambda x: 'iter_{}_'.format(i) in x and x.endswith('smi'), os.listdir(res_pn)), []),
	xrange(1, 8)
)
print(', '.join(map(str, map(len, results))))

print('Scanning strings...')
for i, res in enumerate(results):
	if s in res:
		print(i)

print('Scanning graphs...')
ref = rdkitmol2graph(Chem.MolFromSmiles(s))
mol_results = map(lambda x: map(rdkitmol2graph, map(Chem.MolFromSmiles, x)), results)
for i, res in enumerate(mol_results):
	if any(map(lambda x: is_isomorph(x, ref), res)):
		print(i)
		break
