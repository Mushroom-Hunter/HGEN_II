# This is the project for 2021 HGEN II
# Author: Wanying Zhu
# Email: wanying.zhu.1@vanderbilt.edu

import multiprocessing # For multiprocessing purpose
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import chisquare
# Can also use this one to calculate chi squared p value
# from scipy.stats import chi2

verbose = True # For debugging, print out variables
save_pca_result = True # Save PCA results for my own plotting, but not required in this assignment.
# fn = 'data/HGEN8341_finalproject_data.txt'
# df = pd.read_csv(fn, sep='\t', dtype='str')
fn = 'data/HGEN8341_finalproject_data.txt.gz'
df = pd.read_csv(fn, sep='\t', compression='gzip', dtype='str')

if verbose:
    print('Starting HGEN II final project by Wanying Zhu')
    print('- Original dataset shape:', df.shape)
    print(df.head())

# ------------------ Q1. Allele frequency estimates ------------------
# In the original data file, 0 represents reference allele, other numbers (1, 2 or more) represent alternative allele

# Step 1. Separate 0/1 (alleles) to 0 and 1
# Create a new dataframe (df_allele_separated) with 2 column for each allele, eg. rs116742944_a1 rs116742944_a2
# then separate genotype like 0/1 into the two columns
# Original datafarme df will be needed later, so a new dataframe is necessary
if verbose:
    print('\n\nQ1. Allele frequency estimates')
    print('- Calculate allele frequency of alternative allele:', )

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
    print('- Separate alleles into two columns, new dataframe shape is:', df_allele_separated.shape)
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

# ------------------ Q2. Hardy‑Weinberg Equilibrium ------------------
# Only calculate HWE for variants with MAF>0.05
df_variants_HWE = pd.read_csv('AF.txt', sep='\t')
threshold_maf = 0.05
mask_drop =(df_variants_HWE['alt_allele_frequency']<=threshold_maf) | (df_variants_HWE['alt_allele_frequency']>=1-threshold_maf)
df_variants_HWE.drop(labels=(df_variants_HWE[mask_drop]).index, inplace=True) # Drop where MAF<=0.05

if verbose:
    print('\n\n=======================================\nQ2. HWE')
    print('- Only calculate HWE for SNPs with MAF >', threshold_maf)
    print('- Number of dropped SNPs due to small MAF: ', mask_drop.shape[0])
    print('- Number of SNPs with MAF>0,05: ', df_variants_HWE.shape[0])

# Step 1. Calculate minor allele frequency
df_variants_HWE['MAF'] = df_variants_HWE['alt_allele_frequency'] # Create a column of MAF, set initial values to 'alt_allele_frequency' column
mask_AF_gt_05 = df_variants_HWE['alt_allele_frequency']>=0.5 # mask of SNPs with allele frequency >0.5

# Only this or df_variants_HWE.loc[mask_AF_gt_05, ['MAF']] works, could be a bug of pandas
# df_variants_HWE.loc[mask_AF_gt_05]['MAF'] does not work, even though it also returns a Series for assignment
df_variants_HWE.loc[mask_AF_gt_05, 'MAF'] = 1 - df_variants_HWE.loc[mask_AF_gt_05]['MAF']
# df_variants_HWE.to_csv('AF_and_alt_AF_MAF_gt_005_only.txt', sep='\t', index=False) # For debugging purpose

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

# Calculate expected frequencies of 0, 1 and 2 using df_pair_count_frequency.
# Note that 0 is ref allele, not necessarily major allele, so should not use 1-MAF
df_variants_HWE['0_count_exp'] = (1 - df_variants_HWE['alt_allele_frequency'])**2 * number_of_individuals
df_variants_HWE['1_count_exp'] = 2 * df_variants_HWE['MAF'] * (1 - df_variants_HWE['MAF']) * number_of_individuals
df_variants_HWE['2_count_exp'] = df_variants_HWE['alt_allele_frequency']**2 * number_of_individuals

# Merge df_variants_HWE and df_pair_count_frequency to merge observed and expected counts into the same dataframe
df_merged_exp_obs_freq = df_variants_HWE.merge(df_pair_count_frequency, left_on='rsID', right_index=True)
# df_merged_exp_obs_freq.to_csv('optional_obs_exp_counts_gt_005_only_for_HWE.txt', sep='\t', index=False) # For debugging

if verbose:
    print('- SNPs with MAF>0.5, observed and expected counts:')
    print(df_merged_exp_obs_freq.head())

# Apply chisquare() to calculated p values
# Note: Use ddof=1 to set degree of freedom to 1, otherwise it will be observation-1=3-1=2
#       Second item in the return value of chisquare() is p value
df_merged_exp_obs_freq['HWE_p_vals'] = df_merged_exp_obs_freq.iloc[:, 3:].apply(lambda x: chisquare(f_exp = x[:3], f_obs = x[-3:], ddof=1)[1], axis=1)

# Step 4. Output into HWE.txt
df_merged_exp_obs_freq[['rsID', 'HWE_p_vals']].to_csv('HWE.txt', sep='\t', index=False)

# ------------------ Q3. Linkage Disequilibrium ------------------
# Can also read in AF.txt file to create this dataframe
df_variants_af = pd.DataFrame({'rsID':lst_variant_IDs,'alt_allele_frequency':lst_alt_allele_frequency})
threshold_maf = 0.05
mask_LD_drop =(df_variants_af['alt_allele_frequency']<threshold_maf) | (df_variants_af['alt_allele_frequency']>1-threshold_maf)
df_variants_af.drop(labels=(df_variants_af[mask_LD_drop]).index, inplace=True) # Drop where MAF<0.05
number_of_SNPS_to_use = 100 # Only calculate LD for the first 100 SNPs (4950 pairs)
if verbose:
    print('\n\n=======================================\nQ3. LD')
    print('- Calculate pairwise D, D\' and R2 of the first', number_of_SNPS_to_use, 'SNPs')

# ----------- Helper functions -----------
# Estimate haplotype frequencies of snp1-snp2 with EM algorithm
def estimate_haplotype_frequency(snp1, snp2):

    return 1

# Calculate LD
def get_LD():
    pass
# ----------- End of helper functions -----------

# Calculate LD (D, D' and R2) for the first 100 pairs of SNPs
df_variants_af.reset_index(drop=True, inplace=True) # Reset index since it was modified due to previous drop MAF<0.05
df_variants_af.drop(labels=[i for i in range(number_of_SNPS_to_use, df_variants_af.shape[0])], inplace=True)

''' Not necessary
# Calculate major and minor allele frequency
df_variants_af['major_af'] = df_variants_af['alt_allele_frequency']
df_variants_af.loc[df_variants_af['major_af']<0.5, ['major_af']] = 1 - df_variants_af['major_af']
df_variants_af['minor_af'] = 1 - df_variants_af['major_af']
'''

lst_SNPs_to_calculate_LD = df_variants_af['rsID'] # Actually this is a Series, not a list
ld_D = 0
ld_D_prime = 0
ld_r2 = 0

c = 0
for i in range(number_of_SNPS_to_use):
    snp1 = lst_SNPs_to_calculate_LD[i]
    for j in range(i, number_of_SNPS_to_use):
        snp2 = lst_SNPs_to_calculate_LD[j]

        # Get allele frequencies for the two SNPs
        snp1_major_af = df_variants_af.iloc[i, 2]
        snp1_minor_af = df_variants_af.iloc[i, 3]
        snp2_major_af = df_variants_af.iloc[j, 2]
        snp2_minor_af = df_variants_af.iloc[j, 3]

        # Get genotype counts of the two SNPs (two-locus genotype)
        # Note:
        # 1. There are 3 genotypes for a single SNP: 0 (00), 1(01) and 2(11),
        #    so for two loci there are 3*3=9 possible two-locus genotypes
        # 2. Explanation of naming below: n_00 means SNP1 and SNP2 are both 00. n_12 means SNP1 is 01 and SNP2 is 11
        # 3. n_11 is the only one with ambiguous haplotype. Since both SNPs are 01, we can have two possible haplotypes:
        #    (1) SNP1-SNP2: 0-1 and 1-0; or (2) SNP1-SNP2: 0-0 and 1-1
        # 4. Use EM algorithm to get haplotype frequency
        n_00 = df[(df[snp1]=='0/0') & (df[snp2]=='0/0')].shape[0] # haplotype (SNP1|SNP2): 00|00
        n_01 = df[(df[snp1]=='0/0') & (df[snp2]=='0/1')].shape[0] # haplotype (SNP1|SNP2): 00|01
        n_02 = df[(df[snp1]=='0/0') & (df[snp2]=='1/1')].shape[0] # haplotype (SNP1|SNP2): 01|01
        n_10 = df[(df[snp1]=='0/1') & (df[snp2]=='0/0')].shape[0] # haplotype (SNP1|SNP2): 00|10
        n_11 = df[(df[snp1]=='0/1') & (df[snp2]=='0/1')].shape[0] # haplotype (SNP1|SNP2): 01|10 or 00|11
        n_12 = df[(df[snp1]=='0/1') & (df[snp2]=='1/1')].shape[0] # haplotype (SNP1|SNP2): 01|11
        n_20 = df[(df[snp1]=='1/1') & (df[snp2]=='0/0')].shape[0] # haplotype (SNP1|SNP2): 10|10
        n_21 = df[(df[snp1]=='1/1') & (df[snp2]=='0/1')].shape[0] # haplotype (SNP1|SNP2): 10|11
        n_22 = df[(df[snp1]=='1/1') & (df[snp2]=='1/1')].shape[0] # haplotype (SNP1|SNP2): 11|11

        # ------- EM to estimate haplotype frequencies -------

        # ---- E step ----
        # Initial counts
        # Calculate expected counts of in n_11, ie. when SNP1-SNP2 is 01/01 (hets in both loci)
        # !!!! I am not completely sure about this one, since initially all haplotype frequencies are unknown
        # !!!! For initial counts, assume half of n_11 is 0|1:1|0, the other half is 0|0:1|1
        # !!!! Probably won't matter since we are "guessing" anyway (yes! I tested with other values)
        n_01_10 = n_11 * 0.5 # 01|10
        n_00_11 = n_11 * 0.5 # 00|11

        # haplotype counts of SNP1-SNP2
        # There are 2*2=4 types of haplotype in total (00, 01, 10, 11)
        # number_of_individuals * 2
        haplotype_count_00 = n_00 * 2 + n_01 + n_10 + n_00_11
        haplotype_count_01 = n_01 + n_02 * 2 + n_01_10 + n_12
        haplotype_count_10 = n_10 + n_01_10 + n_20 * 2 + n_21
        haplotype_count_11 = n_00_11 + n_12 + n_21 + n_22 * 2

        if verbose:
            print('\n-', snp1, snp2)
            print(' n_00\tn_01\tn_02\tn_10\tn_11\tn_12\tn_20\tn_21\tn_22')
            print(' '+str(n_00), n_01, n_02, n_10, n_11, n_12, n_20, n_21, n_22, sep='\t\t')
            print(' Initial expected haplotype counts (0|0, 0|1, 1|0, 1|1):',
                  haplotype_count_00, haplotype_count_01, haplotype_count_10, haplotype_count_11)

        # ---- M step ----
        # Update haplotype frequencies with new haplotype counts
        # Stop iteration when changes <10^-5, or reach 1 million iterations
        # From testing, estimated haplotype frequencies usually converge within 20 iterations in this assignment
        for i in range(1000000):
            # Use updated haplotype frequencies to calculate expected counts within n_11
            hap_freq_00 = haplotype_count_00 / (number_of_individuals * 2)
            hap_freq_01 = haplotype_count_01 / (number_of_individuals * 2)
            hap_freq_10 = haplotype_count_10 / (number_of_individuals * 2)
            hap_freq_11 = haplotype_count_11 / (number_of_individuals * 2)
            n_01_10 = n_11 * (2*hap_freq_01*hap_freq_10)/(2*hap_freq_01*hap_freq_10 + 2*hap_freq_00*hap_freq_11) # 01|10
            n_00_11 = n_11 - n_01_10 # 00|11

            print(hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11)

            # Stop iteration when converged
            if i==0:
                # Use prev_hap_freq to track changes of estimated haplotype. Stop iteration when changes are too small
                prev_hap_freq_00 = haplotype_count_00 / (number_of_individuals * 2)
                prev_hap_freq_01 = haplotype_count_01 / (number_of_individuals * 2)
                prev_hap_freq_10 = haplotype_count_10 / (number_of_individuals * 2)
                prev_hap_freq_11 = haplotype_count_11 / (number_of_individuals * 2)
            elif abs(prev_hap_freq_00 - hap_freq_00) < 1e-5 or \
                abs(prev_hap_freq_01 - hap_freq_01) < 1e-5 or \
                abs(prev_hap_freq_10 - hap_freq_10) < 1e-5 or \
                abs(prev_hap_freq_11 - hap_freq_11) < 1e-5:
                break

            haplotype_count_00 = n_00 * 2 + n_01 + n_10 + n_00_11
            haplotype_count_01 = n_01 + n_02 * 2 + n_01_10 + n_12
            haplotype_count_10 = n_10 + n_01_10 + n_20 * 2 + n_21
            haplotype_count_11 = n_00_11 + n_12 + n_21 + n_22 * 2

            prev_hap_freq_00 = hap_freq_00
            prev_hap_freq_01 = hap_freq_01
            prev_hap_freq_10 = hap_freq_10
            prev_hap_freq_11 = hap_freq_11

        c = c+1
        if c>10: break
    if c>10: break




# ------------------ Q4. Principal component analysis ------------------
if verbose:
    print('\n\n=======================================\nQ4. PCA')
# Create a new dataframe (df_additive_genotype) from original dataframe df, convert 0/0, 0/1, 1/1 to 0, 1, 2
df_additive_genotype = df.replace(to_replace='0/0', value=0).copy()
df_additive_genotype.replace(to_replace='0/1', value=1, inplace=True)
df_additive_genotype.replace(to_replace='1/1', value=2, inplace=True)

if verbose:
    print('- Replace 0/0, 0/1, 1/1 with 0, 1, 2 in a new DataFrame (df_additive_genotype):')
    print(df_additive_genotype.head())

# Perform PCA
pca_result = PCA() # Keep the 150 PCs by default (since sample size=150, smaller than SNP size)
pc_vals = pca_result.fit_transform(df_additive_genotype.iloc[:, 2:]) # Fit a PCA model and calculate PCs
if verbose:
    print('- Perform PCA, % variance explained by the first 10 PCs are:')
    for var in pca_result.explained_variance_ratio_[:10]*100:
        print('\t{:.2f}%'.format(var), end='')
    print()

if save_pca_result:
    # Save PC values
    df_PCs = pd.DataFrame(pc_vals)
    df_PCs.columns = ['PC'+str(x+1) for x in df_PCs.columns]
    df_PCs.to_csv('PCs.txt', sep='\t',index=False)

    # Save variance explained by each PC
    with open('PCA_variance_ratio_explained_by_each_PC.txt', 'w') as fh:
        fh.write('PC_number' + '\t' + 'variance_explained' + '\t' + 'variance_explained_ratio' + '\n')
        for i in range(len(pca_result.explained_variance_ratio_)):
            fh.write(str(i+1)+'\t'+str(pca_result.explained_variance_[i])+'\t'+str(pca_result.explained_variance_ratio_[i])+'\n')
