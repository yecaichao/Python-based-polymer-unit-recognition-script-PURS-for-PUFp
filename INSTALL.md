# Installation and Deployment

This repository was validated with the packaged sample file `test.csv`.

## Create the environment

Use the environment file:

```bash
conda env create -f environment_purs.yml
conda activate purs-review-py39
```

## Run the workflow

```bash
python get-polymer-unit.py test.csv --output-dir output
python polymer-unit-classify.py --ring-total-list output/ring_total_list.csv --index-data output/index_data.csv --output-dir output
```

## Input file

The included `test.csv` file can be used as a sample input.
The input csv should contain at least:

- `name`
- `smiles`

## Output files

Running `get-polymer-unit.py` produces:

- `ring_total_list.csv`
- `one_hot.csv`
- `number.csv`
- `adjacent_matrix.csv`
- `node_matrix.csv`
- `index_data.csv`

Running `polymer-unit-classify.py` produces:

- `ring_df.csv`
- `type_frame.csv`

## Notes

- `get-polymer-unit.py` accepts `--output-dir` so generated files can be written into a separate folder.
- `polymer-unit-classify.py` accepts `--ring-total-list`, `--index-data`, and `--output-dir`.
- The workflow keeps the original logic and output format.
