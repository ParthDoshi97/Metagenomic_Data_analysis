## Metagenomic Data Analysis and Visualization
![License](https://img.shields.io/badge/license-MIT-green)
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![R version](https://img.shields.io/badge/r-4.2.1-red.svg)


This repository contains R and Python scripts for performing comprehensive metagenomic data analysis. The scripts are designed to calculate Alpha & Beta Diversity and create Taxonomy pie charts from ASV and OTU data collected from 16s RNA sequences. This script was were create to create specfic data analysis for one of my project 

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
install.packages(c("phyloseq", "vegan", "ape", "ggplot2", "reshape2", "dplyr", "ggpubr"))
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
source("Alpha_&_Beta_diversity.R")
```
Running the Python Script
To perform additional analysis and visualization with the Python script:

```bash
python pie_chart.py
```
## Example plot the script will create

#### Figure 1: Alpha Diversity
<img src="https://github.com/ParthDoshi97/Metagenomic_Data_analysis/assets/126096120/6845e5c5-b757-4213-a409-97d5b91abc1a" alt="Alpha Diversity" width="500" height="400">

#### Figure 2: Beta Diversity
<img src="https://github.com/ParthDoshi97/Metagenomic_Data_analysis/assets/126096120/e978596a-f48a-46bc-aa8a-0a3800102d31" alt="Beta Diversity" width="500" height="400">

#### Figure 3: Taxonomy Pie Chart
<img src="https://github.com/ParthDoshi97/Metagenomic_Data_analysis/assets/126096120/47936fdc-ef12-45ff-b0c2-46dec6a8ff6d" alt="Taxonomy Pie Chart" width="600" height="500">

## License
This project is licensed under the MIT License - see the LICENSE file for detail
