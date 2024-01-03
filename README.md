# Thiophene Reactivity toward Singlet Oxygen
GitHub repository guide for predicting thiophene reactivity with singlet oxygen, covering solvent-free and methanol environments. Instructions include file preparation, ORCA installation, structure optimization, and running the script for results.

# Introduction

Welcome to our GitHub repository! This README provides instructions for users who are interested in executing our script to predict the reactivity of thiophene in the presence of singlet oxygen. The guide is designed to assist you through each step, ensuring a successful execution of the script, whether your focus is on reactivity predictions without solvent effects or within a methanol environment.

## Prerequisites

Before starting, ensure that you have downloaded the code from this repository and have python installed.

## Getting Started

### Prepare Your Inputs:

- Inside the "Input for ORCA" folder, find two sub-folders: "Gas" and "Methanol".
- Choose "Gas" if your objective is to predict thiophene reactivity with singlet oxygen without solvent effect. For predictions in a methanolic environment, use the inputs in "Methanol".

### File Format Requirements:

- ORCA works with `.xyz` files.
- Use Avogadro or another editing and visualization tool to create `.xyz` files.
- If your tool doesn’t support `.xyz` file creation, create a `.mol` file, place it in the "convert" folder, and run `Run_convert.bat` to convert it.

### Charge Consideration:

- Ensure that your designed thiophene compounds always have a net charge of zero.
- Include the counter-ion for ionic compounds.

## Execution Steps

### Install ORCA:

- Install ORCA version 5.0.3 (recommended) for compatibility.

### Optimize Electronic Structure:

- Optimize the electronic structure of thiophene in either the gas phase or methanol environment.

### Calculate Single Points:

- Conduct single-point calculations for the thiophene with N, N+1, and N-1 electron states in the chosen phase (gas or methanol).

### Run the Script:

- Place the outputs in the downloaded folder.
- Choose `Run_ideal(no solvent effect).bat` or `Run_in_Methanol.bat` depending on the environment.

### View Results:

- The reactivity predictions of thiophene with singlet oxygen will be displayed on the screen.

## Special Feature

The script is specifically designed to identify the structure of thiophene and calculate its reactivity with singlet oxygen, even for compounds containing multiple thiophene units.

## Conclusion

Following these steps, you will be able to successfully predict the reactivity of thiophene in the presence of singlet oxygen. For any issues, please raise an issue in the repository for assistance. Happy computing!

