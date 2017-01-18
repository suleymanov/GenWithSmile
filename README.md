# GenWithSmile
SMILES generator for small molecules repertoire

**Install**

pip install -r requirements.txt

**Run**

1) One core

python process.py --one_core sample_inputs/one_core.json

2) Linkers

python process.py --source_target sample_inputs/source_target.json

3) Combinations

python process.py --combinations sample_inputs/combinations.json

See sample_inputs/INFO for configs description and gws/tests/demos for several use cases.
