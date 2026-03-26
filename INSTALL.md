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
python get-polymer-unit.py
python polymer-unit-classify.py
```

## Input file

The included `test.csv` file can be used as a sample input.

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

- The current scripts expect the input csv file to be in the same folder.
- The workflow keeps the original logic and output format.
