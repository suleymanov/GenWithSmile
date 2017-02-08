from gws.src.io.utils import IOUtils


def write_config(fn, config):
	IOUtils.write_json(fn, config)
