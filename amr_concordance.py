import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import cohen_kappa_score, confusion_matrix
from math import pi
import os
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Argument parser for dynamic input
parser = argparse.ArgumentParser(description="Dynamic AMR Concordance Analysis Script")
parser.add_argument("--input", required=True, help="Path to the input AMR table CSV file")
parser.add_argument("--output", default="output_dynamic", help="Directory to save the output plots and metrics")
parser.add_argument("--phenotypic_col_carb", default="Phenotypic Carbapenem susceptibility",
                    help="Column name for phenotypic carbapenem susceptibility")
parser.add_argument("--phenotypic_col_coli", default="Phenotypic Colistin susceptibility",
                    help="Column name for phenotypic colistin susceptibility")
parser.add_argument("--gene_mutation_col", default="Gene Mutation",
                    help="Column name for colistin gene mutations")
parser.add_argument("--carb_genes_col", default="Carbapenem Resistance Genes",
                    help="Column name for carbapenem resistance genes")
parser.add_argument("--exclude_carb_genes", nargs="*", default=["NDM-5", "ompF_Q88STOP"],
                    help="List of carbapenem resistance genes to exclude")
parser.add_argument("--exclude_colistin_mutations", nargs="*", default=["pmrB_Y358N", "pmrB_E123D"],
                    help="List of colistin mutations to exclude")

args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output, exist_ok=True)

# Load the data
logging.info(f"Loading data from {args.input}")
try:
    amr_data = pd.read_csv(args.input)
except FileNotFoundError:
    logging.error(f"File {args.input} not found!")
    exit(1)

# Replace "I" with "S" in phenotypic susceptibility columns
if args.phenotypic_col_carb in amr_data.columns:
    amr_data[args.phenotypic_col_carb] = amr_data[args.phenotypic_col_carb].replace("I", "S")
else:
    logging.error(f"Phenotypic column {args.phenotypic_col_carb} not found in the input data!")
    exit(1)

if args.phenotypic_col_coli in amr_data.columns:
    amr_data[args.phenotypic_col_coli] = amr_data[args.phenotypic_col_coli].replace("I", "S")
else:
    logging.error(f"Phenotypic column {args.phenotypic_col_coli} not found in the input data!")
    exit(1)

# Classification functions
def classify_carbapenem(row, exclude_genes):
    phenotypic_resistance = "R" if row[args.phenotypic_col_carb] == "R" else "S"
    genes = str(row[args.carb_genes_col]).split("&") if pd.notnull(row[args.carb_genes_col]) else []
    genotypic_resistance = "R" if any(gene not in exclude_genes for gene in genes) else "S"
    return phenotypic_resistance, genotypic_resistance

def classify_colistin(row, exclude_substitutions):
    phenotypic_resistance = "R" if row[args.phenotypic_col_coli] == "R" else "S"
    gene_mutation = str(row[args.gene_mutation_col]) if pd.notnull(row[args.gene_mutation_col]) else "-"
    if exclude_substitutions and gene_mutation in args.exclude_colistin_mutations:
        genotypic_resistance = "S"  # Mark these as non-resistant
    else:
        genotypic_resistance = "R" if gene_mutation != "-" else "S"
    return phenotypic_resistance, genotypic_resistance

# Analyze data
results_dict = {
    "phenotypes_carbapenem": [], "genotypes_carbapenem": [],
    "phenotypes_colistin_incl": [], "genotypes_colistin_incl": [],
    "phenotypes_colistin_excl": [], "genotypes_colistin_excl": []
}

for _, row in amr_data.iterrows():
    carb_pheno, carb_geno = classify_carbapenem(row, args.exclude_carb_genes)
    results_dict["phenotypes_carbapenem"].append(carb_pheno)
    results_dict["genotypes_carbapenem"].append(carb_geno)
    
    coli_pheno_incl, coli_geno_incl = classify_colistin(row, exclude_substitutions=False)
    results_dict["phenotypes_colistin_incl"].append(coli_pheno_incl)
    results_dict["genotypes_colistin_incl"].append(coli_geno_incl)
    
    coli_pheno_excl, coli_geno_excl = classify_colistin(row, exclude_substitutions=True)
    results_dict["phenotypes_colistin_excl"].append(coli_pheno_excl)
    results_dict["genotypes_colistin_excl"].append(coli_geno_excl)

# Calculate metrics
def calculate_metrics(phenotypes, genotypes):
    tn, fp, fn, tp = confusion_matrix(phenotypes, genotypes, labels=["S", "R"]).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    return sensitivity, specificity, ppv, npv

metrics = {}
for key in ["carbapenem", "colistin_incl", "colistin_excl"]:
    pheno_key = f"phenotypes_{key}"
    geno_key = f"genotypes_{key}"
    kappa = cohen_kappa_score(results_dict[pheno_key], results_dict[geno_key])
    sens, spec, ppv, npv = calculate_metrics(results_dict[pheno_key], results_dict[geno_key])
    metrics[key] = [kappa, sens, spec, ppv, npv]

# Save metrics
metrics_file = os.path.join(args.output, "concordance_metrics_dynamic.txt")
with open(metrics_file, "w") as f:
    for key, values in metrics.items():
        f.write(f"{key.capitalize()} Metrics:\n")
        for metric_name, value in zip(["Cohen's Kappa", "Sensitivity", "Specificity", "PPV", "NPV"], values):
            f.write(f"{metric_name}: {value:.4f}\n")
        f.write("\n")
logging.info(f"Metrics saved to {metrics_file}")

# Bar Plot
categories = ["Cohen's Kappa", "Sensitivity", "Specificity", "PPV", "NPV"]
carbapenem_values = metrics["carbapenem"]
colistin_incl_values = metrics["colistin_incl"]
colistin_excl_values = metrics["colistin_excl"]

x = np.arange(len(categories))
bar_width = 0.25

fig, ax = plt.subplots(figsize=(12, 7))
bars1 = ax.bar(x - bar_width, carbapenem_values, bar_width, label='Carbapenem')
bars2 = ax.bar(x, colistin_incl_values, bar_width, label='Colistin (incl. mutations)')
bars3 = ax.bar(x + bar_width, colistin_excl_values, bar_width, label='Colistin (excl. mutations)')

ax.set_ylabel('Values')
ax.set_title('Concordance Metrics Comparison')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()

for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')

plt.tight_layout()
bar_plot_path = os.path.join(args.output, "concordance_metrics_comparison_bar_dynamic.pdf")
plt.savefig(bar_plot_path)
plt.close()
logging.info(f"Bar plot saved to {bar_plot_path}")

# Radar Plot
angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
angles += angles[:1]

fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(polar=True))
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)
plt.xticks(angles[:-1], categories)
plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], ["0.2", "0.4", "0.6", "0.8", "1"], color="grey", size=7)
plt.ylim(0, 1)

# Plot each dataset
def add_radar_plot(values, label, color):
    values += values[:1]  # Close the radar chart
    ax.plot(angles, values, linewidth=2, linestyle='solid', label=label, color=color)
    ax.fill(angles, values, color=color, alpha=0.25)

add_radar_plot(carbapenem_values, "Carbapenem", "blue")
add_radar_plot(colistin_incl_values, "Colistin (incl. mutations)", "orange")
add_radar_plot(colistin_excl_values, "Colistin (excl. mutations)", "green")

plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
radar_plot_path = os.path.join(args.output, "concordance_metrics_comparison_radar_dynamic.pdf")
plt.savefig(radar_plot_path)
plt.close()
logging.info(f"Radar plot saved to {radar_plot_path}")
