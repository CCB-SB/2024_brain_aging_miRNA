#!/bin/bash

# preprocessing
snakemake -s preprocessing.Snakefile --use-conda --conda-frontend conda --configfile configs/config_preprocessing_CA1.yaml -j20 -k # works
snakemake -s preprocessing.Snakefile --use-conda --conda-frontend conda --configfile configs/config_preprocessing_CA2_diet_restriction.yaml -j20 -k # works
snakemake -s preprocessing.Snakefile --use-conda --conda-frontend conda --configfile configs/config_preprocessing_CA2_plasma_injection.yaml -j20 -k #works
snakemake -s preprocessing.Snakefile --use-conda --conda-frontend conda --configfile configs/config_preprocessing_microglia.yaml -j20 -k # works

python scripts/isomir_preparation.py "raw_data/data_CA1_aging/isomiR_expression_raw_aging.tsv" "raw_data/data_CA1_aging/miRNA_filtered_quantification_rpmm_norm.detection_rate_10p_per_brain_region.csv" "data/mapping_infos/mapping_vs_mm10_mis1_CA1_CA2.csv" "data/annotations/annotation_CA1_aging_filtered.csv" "raw_data/data_CA1_aging/iso_filtered_quantification_rpmm_norm.detection_rate_10p_per_brain_region.csv" #works

# python create some expressions for only male or female or young
python scripts/select_samples_from_table.py "CA1_aging_age=3m_12m_15m" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_age=3m_12m_15m/" "data/annotations" "miRNA" # works
python scripts/select_samples_from_table.py "CA1_aging_age=3m_12m_15m_male" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_age=3m_12m_15m_male/" "data/annotations" "miRNA" # works
python scripts/select_samples_from_table.py "CA1_aging_age=3m_12m_15m_female" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_age=3m_12m_15m_female/" "data/annotations" "miRNA" # works

python scripts/select_samples_from_table.py "CA1_aging_male" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_male/" "data/annotations" "miRNA" "tRNA" # works
python scripts/select_samples_from_table.py "CA1_aging_without_age=26m_28m_female" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_without_age=26m_28m_female/" "data/annotations" "miRNA" "tRNA" # works

python scripts/select_samples_from_table.py "CA1_aging_without_age=26m_28m" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_without_age=26m_28m/" "data/annotations" "miRNA" # works

python scripts/select_samples_from_table.py "CA1_aging_without_pon_male" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_without_pon_male/" "data/annotations" "miRNA" # works
python scripts/select_samples_from_table.py "CA1_aging_without_age=26m_28m_without_pon_female" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_without_age=26m_28m_without_pon_female/" "data/annotations" "miRNA" # works

python scripts/select_samples_from_table.py "CA1_aging_without_age=26m_28m_without_pon" "raw_data/data_CA1_aging/" "raw_data/data_CA1_aging_without_age=26m_28m_without_pon/" "data/annotations" "miRNA" # works

# aging cohorte all RNA classes
snakemake -s rna_classes.Snakefile --use-conda --conda-frontend conda --configfile configs/config_rna_classes_CA1.yaml -j20 -k # works

# aging cohorte tRNA
snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_tRNA_CA1.yaml -j20 -k # works
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_tRNA_CA1_sex_specific_male.yaml -j20 -k # works
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_tRNA_CA1_sex_specific_female.yaml -j20 -k #works
snakemake -s sex_specific_male_female.Snakefile --use-conda --conda-frontend conda --configfile configs/config_tRNA_CA1_sex_specific_male_female.yaml -j20 -k # works

# aging cohorte miRNA
snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_young.yaml -j20 -k #works  
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_young_male.yaml -j20 -k #works
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_young_female.yaml -j20 -k #works
snakemake -s sex_specific_male_female.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_young_male_female.yaml -j20 -k #works
snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1.yaml -j20 -k # works
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_male.yaml -j20 -k #works
snakemake -s sex_specific.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_female.yaml -j20 -k #works
snakemake -s sex_specific_male_female.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_male_female.yaml -j20 -k #works
snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA1_without_old.yaml -j20 -k #works
snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_isomiR_CA1.yaml -j20 -k #works

# diet restriction experiment miRNA
snakemake -s rna_classes.Snakefile --use-conda --conda-frontend conda --configfile configs/config_rna_classes_CA2_diet_restriction.yaml -j20 -k #works
snakemake -s CA2.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA2_diet_restriction.yaml -j20 -k #works

# young plasma in old mice injection experiment miRNA
snakemake -s rna_classes.Snakefile --use-conda --conda-frontend conda --configfile configs/config_rna_classes_CA2_plasma_injection.yaml -j20 -k #works
snakemake -s CA2.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_CA2_plasma_injection.yaml -j20 -k #works

# microglia miRNA
snakemake -s rna_classes.Snakefile --use-conda --conda-frontend conda --configfile configs/config_rna_classes_microglia.yaml -j20 -k #works
snakemake -s microglia.Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_microglia.yaml -j20 -k #works

# Opional: mRNA miRNA interplay (detail from Hahn et al.)
# Check if the mrna variable is set
if [ -z "$MRNA" ]; then
    echo "Error: MRNA variable is not set. Please set it to 1 or 0."
    exit 1
fi

# Execute commands based on the value of mrna
if [ "$MRNA" -eq 1 ]; then
    echo "mRNA is set to 1. Executing commands for MRNA = 1..."
    # Add commands to execute when MRNA=1
    echo "Running analysis for MRNA"
    snakemake -s Snakefile --use-conda --conda-frontend conda --configfile configs/config_miRNA_mRNA_interplay_CA1.yaml -j20 -k #works

elif [ "$MRNA" -eq 0 ]; then
    echo "MRNA is set to 0. Executing commands for MRNA = 0..."
    # Add commands to execute when MRNA=0
    echo "Skipping MRNA analysis"
else
    echo "Error: Invalid value for MRNA. Please set it to either 1 or 0."
    exit 1
fi

#Opional: ROSMAP human brain miRNA
# Check if the mrna variable is set
if [ -z "$HUMAN" ]; then
    echo "Error: HUMAN variable is not set. Please set it to 1 or 0."
    exit 1
fi

# Execute commands based on the value of HUMAN
if [ "$HUMAN" -eq 1 ]; then
    echo "HUMAN is set to 1. Executing commands for HUMAN = 1..."
    # Add commands to execute when HUMAN=1
    echo "Running analysis for HUMAN"
    snakemake -s preprocessing.Snakefile --use-conda --conda-frontend conda --configfile configs/config_preprocessing_human.yaml -j20 -k # works
    snakemake -s human.Snakefile --use-conda --conda-frontend conda --configfile configs/config_human.yaml -j20 -k #works

elif [ "$HUMAN" -eq 0 ]; then
    echo "HUMAN is set to 0. Executing commands for HUMAN = 0..."
    # Add commands to execute when HUMAN=0
    echo "Skipping HUMAN analysis"
else
    echo "Error: Invalid value for HUMAN. Please set it to either 1 or 0."
    exit 1
fi

python scripts/extract_paper_figures.py
