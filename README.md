# PURS for PUFp v1.1

This repository provides a lightweight workflow for polymer-unit recognition
and polymer-unit fingerprint (PUFp) generation.

## Included files

- `get-polymer-unit.py`
- `polymer-unit-classify.py`
- `structure_identity_tool.py`
- `test.csv`

## Sample input

`test.csv` is included as a runnable sample input file.

## Quick start

Place the scripts and the input csv in the same folder, then run:

```bash
python get-polymer-unit.py
python polymer-unit-classify.py
```

## Outputs from `get-polymer-unit.py`

- `ring_total_list.csv`
- `one_hot.csv`
- `number.csv`
- `adjacent_matrix.csv`
- `node_matrix.csv`
- `index_data.csv`

## Outputs from `polymer-unit-classify.py`

- `ring_df.csv`
- `type_frame.csv`

## Installation

For a validated environment and deployment instructions, see:

- `INSTALL.md`
- `environment_purs.yml`

## Notes

- Keep the input csv in the same folder as the scripts unless you modify the file paths in the scripts.
- The current scripts preserve the original workflow logic and output structure.

## Citation

If you use PURS in your research, please cite:

[1] Xinyue Zhang, Ye Sheng, Xiumin Liu, Jiong Yang, William A. Goddard III, Caichao Ye*, Wenqing Zhang*. Polymer-unit Graph: Advancing Interpretability in Graph Neural Network Machine Learning for Organic Polymer Semiconductor Materials. J. Chem. Theory Comput., 2024, 20(7), 2908-2920.  
[2] Xinyue Zhang, Genwang Wei, Ye Sheng, Wenjun Bai, Jiong Yang, Wenqing Zhang*, Caichao Ye*. Polymer-Unit Fingerprint (PUFp): An Accessible Expression of Polymer Organic Semiconductors for Machine Learning. ACS Appl. Mater. Interfaces, 2023, 15(17), 21537-21548.  
[3] Caichao Ye, Tao Feng, Weishu Liu*, Wenqing Zhang*. Functional Unit: A New Perspective on Materials Science Research Paradigms. Acc. Mater. Res., 2025, 6(8), 914-920.
