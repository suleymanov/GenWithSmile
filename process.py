import sys

from gws import (
	OneCoreClient,
	SourceTargetClient,
	CombinationsClient,
	BlocksClient,
	SynthpathClient
)


option_matcher = {
	'--one_core': OneCoreClient,
	'--source_target': SourceTargetClient,
	'--combinations': CombinationsClient,
	'--blocks': BlocksClient,
	'--synthpath': SynthpathClient
}


def _process_concept(concept_client, config_fn):
	client = concept_client.from_file(config_fn)
	client.process()


def main(args):
	"""
	Main entry point.
	:param args: list of str
	:return: None
	"""
	assert len(args) == 2

	option = args[0]
	config_fn = args[1]

	concept_client = option_matcher.get(option, None)
	if concept_client:
		_process_concept(concept_client, config_fn)
	else:
		raise KeyError('Unknown option.')


if __name__ == '__main__':
	main(sys.argv[1:])
