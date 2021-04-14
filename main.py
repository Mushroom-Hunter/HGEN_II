# This is the project for 2021 HGEN II
# Author: Wanying Zhu
# Email: wanying.zhu.1@vanderbilt.edu

import pandas as pd
from scipy.stats import chisquare
# Can also use this one to calculate chi squared p value
# from scipy.stats import chi2

verbose = True # For debugging, print out variables
fn = 'data/HGEN8341_finalproject_data.txt'
df = pd.read_csv(fn, sep='\t')

if verbose:
    print('Starting HGEN II final project by Wanying Zhu')
    print('- Original dataset shape:', df.shape)
    print(df.head())

# -------------- Q1. Allele frequency estimates --------------
# In the original data file, 0 represents reference allele, other numbers (1, 2 or more) represent alternative allele

# Step 1. Separate 0/1 (alleles) to 0 and 1
# Create a new dataframe (df_allele_separated) with 2 column for each allele, eg. rs116742944_a1 rs116742944_a2
# then separate genotype like 0/1 into the two columns
# Original datafarme df will be needed later, so a new dataframe is necessary
if verbose:
    print('\n\nQ1. Allele frequency estimates')
    print('Calculate allele frequency of the alternative allele:', )
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
    print('- Separate alleles into tow columns, new dataframe shape is:', df_allele_separated.shape)
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
# Plotting is in a separated file plotting.py
with open('AF.txt', 'w') as fh_out:
    fh_out.write('rsID\talt_allele_frequency\n')
    lst_variant_IDs = df.columns[2:] # Skip ID and status in the header
    for i in range(len(lst_variant_IDs)):
        line = lst_variant_IDs[i] + '\t' + str(lst_alt_allele_frequency[i]) + '\n'
        fh_out.write(line)

# -------------- Q2. Hardy‑Weinberg Equilibrium --------------
# Only calculate HWE for variants with MAF>0.05
df_variants_HWE = pd.read_csv('AF.txt', sep='\t')
threshold_maf = 0.05
mask_drop =(df_variants_HWE['alt_allele_frequency']<=threshold_maf) | (df_variants_HWE['alt_allele_frequency']>=1-threshold_maf)
df_variants_HWE.drop(labels=(df_variants_HWE[mask_drop]).index, inplace=True) # Drop where MAF<=0.05

if verbose:
    print('\n\n=======================================\nQ2. HWE')
    print('- Only calculate HWE for SNPs with MAF >', threshold_maf)
    print('- Dropped number of SNPs due to small MAF: ', mask_drop.shape[0])
    print('- Number of SNPs with MAF>0,05: ', df_variants_HWE.shape[0])

# Step 1. Calculate minor allele frequency
df_variants_HWE['MAF'] = df_variants_HWE['alt_allele_frequency'] # Create a column of MAF, set initial values to 'alt_allele_frequency' column
mask_AF_gt_05 = df_variants_HWE['alt_allele_frequency']>=0.5 # mask of SNPs with allele frequency >0.5
df_variants_HWE.loc[mask_AF_gt_05, 'MAF'] = 1 - df_variants_HWE.loc[mask_AF_gt_05]['MAF'] # First half is a DataFrame, second is a Series, only this way works

# Step 2. Count and calculate frequencies of 0/0, 0/1 (could be 1/0) and 1/1 in original DataFrame df
mask_00 = df.iloc[:,2:]=='0/0' # Skip the first two columns (ID and status)
mask_01 = (df.iloc[:,2:]=='0/1') | (df.iloc[:,2:]=='1/0')
mask_11 = df.iloc[:,2:]=='1/1'

# Store observed frequencies of 0/0 (0), 0/1 (1) and 1/1 (2) to this dataframe
df_pair_count_frequency = pd.concat([mask_00.sum(), mask_01.sum(), mask_11.sum()], axis=1)
df_pair_count_frequency.columns = ['0_count_obs', '1_count_obs', '2_count_obs']
'''
# Do not need to calculate frequencies
df_pair_count_frequency['0_freq_obs'] = df_pair_count_frequency['0_count_obs']/number_of_individuals # divided by number of subjects
df_pair_count_frequency['1_freq_obs'] = df_pair_count_frequency['1_count_obs']/number_of_individuals
df_pair_count_frequency['2_freq_obs'] = df_pair_count_frequency['2_count_obs']/number_of_individuals
'''

if verbose:
    print('- Count and frequency of 0 (0/0), 1 (0/1) and 2 (1/1) for all variants:')
    print(df_pair_count_frequency.head())

# Step 3. Calculate HWE p values
# Note: since minor and major allele frequencies follow HWE,
#       we use 1 as degree of freedom even though there are 3 values 0/0, 0/1 and 1/1
# Use df_variants_HWE and df_pair_count_frequency in this step
# - df_variants_HWE: contains MAFs for SNPs with MAF >0.05
# - df_pair_count_frequency: contains observed counts and frequencies of 0/0, 0/1 and 1/1

# Calculate expected frequencies of 0, 1 and 2 using df_pair_count_frequency
df_variants_HWE['0_count_exp'] = (1 - df_variants_HWE['MAF'])**2 * number_of_individuals
df_variants_HWE['1_count_exp'] = 2 * df_variants_HWE['MAF'] * (1 - df_variants_HWE['MAF']) * number_of_individuals
df_variants_HWE['2_count_exp'] = df_variants_HWE['MAF']**2 * number_of_individuals

# Merge df_variants_HWE and df_pair_count_frequency
df_merged_exp_obs_freq = df_variants_HWE.merge(df_pair_count_frequency, left_on='rsID', right_index=True)

if verbose:
    print('- SNPs with MAF>0.5, observed and expected counts:')
    print(df_merged_exp_obs_freq.head())

# Apply chisquare() to calculated p values
# Note: Use ddof=1 to set degree of freedom to 1, otherwise it will be observation-1=3-1=2
#       Second item in the return value of chisquare() is p value
df_merged_exp_obs_freq['HWE_p_vals']= df_merged_exp_obs_freq.iloc[:, 3:].apply(lambda x: chisquare(f_exp = x[:3], f_obs = x[-3:], ddof=1)[1], axis=1)

# Step 4. Output into HWE.txt
df_merged_exp_obs_freq[['rsID', 'HWE_p_vals']].to_csv('HWE.txt', sep='\t', index=False)

# -------------- Q3. Linkage Disequilibrium --------------
if verbose:
    print('\n\n=======================================\nQ3. LD')



# -------------- Q4. Principal component analysis --------------
if verbose:
    print('\n\n=======================================\nQ4. PCA')