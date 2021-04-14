import matplotlib.pyplot as plt
import qq_plot # This is a in-house code to make QQ plot

# -------------- Question 1. Plot a histogram of allele frequency --------------
# Read in AF.txt file
fn_af = 'AF.txt'
lst_alt_allele_frequency = []
with open(fn_af) as fh:
    line = fh.readline() # Skip header line
    line = fh.readline().strip()
    while line != '':
        lst_alt_allele_frequency.append(float(line.split()[1]))
        line = fh.readline().strip()

fig_af, ax_af = plt.subplots(figsize=(10,5), dpi=150)
ax_af.hist(lst_alt_allele_frequency, bins=50, rwidth=0.8)
ax_af.set_title('Distribution of estimated AF of alternative allele')
ax_af.set_xlabel('Allele frequency')
ax_af.set_ylabel('Count')
fig_af.savefig('Q1_AF_histogram.jpeg')

# Set threshold of allele frequency to plot
threshold_af = 0.05
lst_alt_allele_frequency_non_zeros = []
for i in lst_alt_allele_frequency:
    if i>threshold_af: lst_alt_allele_frequency_non_zeros.append(i)
fig_af_non_zeros, ax_af_non_zeros = plt.subplots(figsize=(10,5), dpi=150)
ax_af_non_zeros.hist(lst_alt_allele_frequency_non_zeros, bins=50, rwidth=0.8)
ax_af_non_zeros.set_title('Distribution of estimated AF of alternative allele (AF>'+str(threshold_af)+')')
ax_af_non_zeros.set_xlabel('Allele frequency')
ax_af_non_zeros.set_ylabel('Count')
fig_af_non_zeros.savefig('Q1_AF_gt_'+ str(threshold_af) +'_histogram.jpeg')


# -------------- Question 2. HWE --------------
# Plot QQ plot for HWE p values
qq_plot.qqplot(filename='HWE.txt', output='Q2_HWE_p_values_QQ_plot.jpeg',
               p_value_column_title = 'HWE_p_vals', title='HWE p values')

# -------------- Question 3. LD --------------



# -------------- Question 4. PCA --------------

