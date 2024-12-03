# AMR Concordance Analysis

This repository provides a dynamic script for analyzing concordance metrics in antimicrobial resistance (AMR) data. The script supports the analysis of phenotypic and genotypic concordance for Carbapenem and Colistin resistance. 

---

## Features
- Dynamically specify input columns for phenotypic and genotypic resistance analysis.
- Exclude specific genes or mutations from Carbapenem or Colistin resistance evaluation.
- Generate publication-ready plots (Bar and Radar charts) and a metrics summary.

---

## Usage

### Prerequisites
- Python 3.7 or above
- Required Python packages:
  - `pandas`
  - `matplotlib`
  - `numpy`
  - `scikit-learn`

Install the required packages using:
```bash
pip install -r requirements.txt
```

### Running the Script
### Example command:
```bash
python amr_concordance.py --input input/AMR_table_demo.csv --output output --exclude_colistin_mutations pmrB_Y358N pmrB_E123D
``` 
--exclude_colistin_mutations pmrB_Y358N pmrB_E123D
--input: Path to the input CSV file containing AMR data.
--output: Directory to save the output plots and metrics.
--exclude_colistin_mutations: List of specific Colistin mutations to exclude.

Additional arguments can be explored in the script using:
```bash
python amr_concordance.py --help
``` 
### Example Input
AMR_table_demo.csv

### Example Outputs

Metrics Summary: concordance_metrics_dynamic.txt
Bar Plot: concordance_metrics_comparison_bar_dynamic.pdf
Radar Plot: concordance_metrics_comparison_radar_dynamic.pdf

---

## Citation
If you are using the amr_concordance.py script, please cite it as follows:

Sharma, V. (2024). amr_concordance.py [Python script]. Retrieved from [https://github.com/vsmicrogenomics/amr_concordance_analysis]


## Acknowledgements
This script utilizes AMRFinderPlus output for extraction. Please acknowledge the use of AMRFinderPlus by referring to the tool at [https://github.com/ncbi/amr]
