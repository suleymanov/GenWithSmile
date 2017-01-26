import os
import sys


def _read_smi(fn):
	with open(fn) as f:
		return f.read().splitlines()


def _show_results_details(args):
	results_pn = args[0]
	# fns = map(
	# 	lambda x: results_pn + os.sep + x, 
	# 	filter(lambda x: x.endswith('smi'), os.listdir(results_pn)))
	# rel_times = map(os.path.getmtime, sorted(fns, key=lambda x: os.path.getmtime(x)))
	# rel_times = map(lambda x: x - rel_times[0], rel_times)
	# print(', '.join(map(str, rel_times)))
	# print(', '.join(map(str, sorted(map(len, map(_read_smi, fns))))))
	fns = sorted(map(
		lambda x: results_pn + os.sep + x,
		filter(lambda x: x.endswith('smi'), os.listdir(results_pn))),
		key=lambda x: os.path.getmtime(x))
	rel_times = map(lambda x: os.path.getmtime(x) - os.path.getmtime(fns[0]), fns)
	# print('\n'.join(map(lambda x: '\t'.join(map(str, x)), zip(fns, rel_times))))
	print('\n'.join(map(lambda x: '\t'.join(map(str, x)), zip(fns, rel_times, map(len, map(_read_smi, fns))))))


def main(args):
	option = args[0]
	if option == '--details':
		return _show_results_details(args[1:])


if __name__ == '__main__':
	main(sys.argv[1:])
