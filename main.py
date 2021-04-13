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
mask_alt_allele = df_allele_separated!='0'
alternative_allele_counts = mask_alt_allele.iloc[:,2:].apply(sum)
number_of_individuals = df.shape[0]
lst_alt_allele_frequency = [] # To store calculated allele frequencies
for i in range(0, len(alternative_allele_counts), 2):
    # Add up counts of two alleles together to get total counts
    total_counts = alternative_allele_counts[i] + alternative_allele_counts[i+1]
    lst_alt_allele_frequency.append(total_counts/(number_of_individuals*2))

# Step 3. Output allele frequencies into 'AF.txt' file
with open('AF.txt', 'w') as fh_out:
    fh_out.write('rsID\talt_allele_frequency')
    lst_variant_IDs = df.columns[2:] # Skip ID and status in the header
    for i in range(len(lst_variant_IDs)):
        line = lst_variant_IDs[i] + '\t' + str(lst_alt_allele_frequency[i]) + '\n'
        fh_out.write(line)

# Step 4. Plot a histogram of allele frequency
fig_af, ax_af = plt.subplots(figsize=(10,5), dpi=150)
ax_af.hist(lst_alt_allele_frequency, bins=50, rwidth=0.8)
ax_af.set_title('Distribution of estimated AF of alternative allele')
ax_af.set_xlabel('Allele frequency')
ax_af.set_ylabel('Count')
fig_af.savefig('Q1_AF_histogram.jpeg')

# If skip zeros
lst_alt_allele_frequency_non_zeros = []
for i in lst_alt_allele_frequency:
    if i>0: lst_alt_allele_frequency_non_zeros.append(i)
fig_af_non_zeros, ax_af_non_zeros = plt.subplots(figsize=(10,5), dpi=150)
ax_af_non_zeros.hist(lst_alt_allele_frequency_non_zeros, bins=50, rwidth=0.8)
ax_af_non_zeros.set_title('Distribution of estimated AF of alternative allele (Non zero AF only)')
ax_af_non_zeros.set_xlabel('Allele frequency')
ax_af_non_zeros.set_ylabel('Count')
fig_af_non_zeros.savefig('Q1_AF_non_zeros_histogram.jpeg')

# -------------- Q2. Hardy‑Weinberg Equilibrium --------------



# -------------- Q3. Linkage Disequilibrium --------------




# -------------- Q4. Principal component analysis --------------