## Metagenomic Data Analysis and Visualization
![License](https://img.shields.io/badge/license-MIT-green)
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![R version](https://img.shields.io/badge/r-4.2.1-red.svg)


This repository contains R and Python scripts for performing comprehensive metagenomic data analysis. The scripts are designed to calculate Alpha & Beta Diversity and create Taxonomy pie charts from ASV and OTU data collected from 16s RNA sequences.

## Features
Alpha Diversity Calculation: Measure the diversity within a single sample.
Beta Diversity Calculation: Measure the diversity between different samples.
Taxonomy Pie Charts: Visualize the composition of microbial communities.

## Installation
To run the scripts, you'll need to have R and Python installed on your system along with the necessary libraries. Follow the instructions below to set up your environment:

### R Setup
Install R from CRAN.
Install the required R packages:

```bash
install.packages(c("vegan", "phyloseq", "ggplot2"))
```
### Python Setup
Install Python from Python.org.
Install the required Python packages using pip:
```bash
pip install numpy pandas matplotlib seaborn
```
## Usage
Running the R Script
To calculate Alpha & Beta Diversity and create Taxonomy pie charts using the R script:

```bash r
source("metagenomic_analysis.R")
```
Running the Python Script
To perform additional analysis and visualization with the Python script:

```bash
python metagenomic_analysis.py
```
## Example plot the script will create

Figure 1: Alpha Diversity


Figure 2: Beta Diversity


Figure 3: Taxonomy Pie Chart


## License
This project is licensed under the MIT License - see the LICENSE file for detail
