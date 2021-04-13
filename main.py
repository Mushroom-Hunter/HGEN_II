# This is the project for 2021 HGEN II
# Author: Wanying Zhu
# Email: wanying.zhu.1@vanderbilt.edu

import pandas as pd
import matplotlib.pyplot as plt

verbose = True # For debugging, print out variables
fn = 'data/HGEN8341_finalproject_data.txt'
df = pd.read_csv(fn, sep='\t')

if verbose:
    print('- Original dataset shape:', df.shape)
    print(df.head())

# -------------- Q1. Allele frequency estimates --------------
# In the original data file, 0 represents reference allele, other numbers (1, 2 or more) represent alternative allele

# Step 1. Separate 0/1 (alleles) to 0 and 1
# Create a new dataframe (df_allele_separated) with 2 column for each allele, eg. rs116742944_a1 rs116742944_a2
# then separate genotype like 0/1 into the two columns
# Original datafarme df will be needed later, so a new dataframe is necessary
if verbose:
    print('\nCalculate allele frequency of the alternative allele:', )
    print(df.head())
lst_separated_column_names = ['ID','status'] # Create column names for new dataframe
for col in df.columns[2:]: # The first two items are ID and status, already included when creating the list
    lst_separated_column_names.append(col+'_a1')
    lst_separated_column_names.append(col + '_a2')

df_allele_separated = pd.DataFrame(columns=lst_separated_column_names)
# Fill ID and status columns
df_allele_separated[lst_separated_column_names[0]] = df[lst_separated_column_names[0]]
df_allele_separated[lst_separated_column_names[1]] = df[lst_separated_column_names[1]]
for i in range(2, len(lst_separated_column_names), 2): # Skip ID and status columns
    df_allele_separated[[lst_separated_column_names[i], lst_separated_column_names[i+1]]] = df[lst_separated_column_names[i]\
                                                                                            .split('_a1')[0]].str.split('/', expand=True)

if verbose:
    print('- Separate alleles into tow columns, new dataset shape is:', df_allele_separated.shape)
    print(df_allele_separated.head())

# Step 2. Calculate alternative allele (ie note '0' alleles) frequency

# -------------- Q2. Hardy‑Weinberg Equilibrium --------------



# -------------- Q3. Linkage Disequilibrium --------------




# -------------- Q4. Principal component analysis --------------