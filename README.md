## HGEN 8341 final project
1. This project is written in python 3  
2. Required modules are **multiprocessing**, **time**, **pandas**, **scipy**, **sklearn**, **matplotlib** and **seaborn**. They should be included if you installed python with Anaconda  
3. For a quick start you can use below statements in command line to generate outputs in the working directory:
	```shell
	python main.py	# Generate text output  
	python plotting.py	# Generate plots
	```
4. Output files are:  
(1) From main.py:  
	```
	# Required files:
	  AF.txt
	  HWE.txt
	  LD_D.txt
	  LD_Dprime.txt
	  LD_r2.txt
	# Optional files (but needed for plotting):
	  PCA_variance_ratio_explained_by_each_PC.txt
	  PCs.txt
	```  
	
(2) From plotting.py:  
	```
	Q1_AF_histogram.jpeg
	Q1_AF_gt_0.05_histogram.jpeg
	Q1_AF_lt_0.01_histogram.jpeg
	Q2_HWE_p_values_QQ_plot.jpeg
	Q2_HWE_p_values_distribution.jpeg
	Q3_LD.jpeg
	Q3_LD_values_and_AF.jpeg
	Q3_LD_values_pairwise_comparison.jpeg
	Q4_PCA.jpeg
	Q4_PCA_colored.jpeg
	```

5. Console output of main.py (on macOS Big Sur with 2-core CPU, finished in 4.91 minutes)
	```console
	# Starting HGEN II final project by Wanying Zhu

	- Original dataset shape: (150, 28746)
		ID status rs555599246  ... rs6140157 rs185180844 rs527420293
	0  NA06984   case         0/0  ...       1/1         0/0         0/0
	1  NA06986   case         0/0  ...       1/1         0/0         0/0
	2  NA06989   case         0/0  ...       1/1         0/0         0/0
	3  NA07000   case         0/0  ...       1/1         0/0         0/0
	4  NA07037   case         0/0  ...       1/1         0/0         0/0

	[5 rows x 28746 columns]


	=======================================
	Q1. Allele frequency estimates
	- Separate alleles into two columns, new dataframe shape is: (150, 57490)
		ID status rs555599246_a1  ... rs185180844_a2 rs527420293_a1 rs527420293_a2
	0  NA06984   case              0  ...              0              0              0
	1  NA06986   case              0  ...              0              0              0
	2  NA06989   case              0  ...              0              0              0
	3  NA07000   case              0  ...              0              0              0
	4  NA07037   case              0  ...              0              0              0

	[5 rows x 57490 columns]
	- Calculate allele frequency of alternative allele


	=======================================
	Q2. HWE
	- Only calculate HWE for SNPs with MAF > 0.05
	- Number of dropped SNPs due to small MAF:  28744
	- Number of SNPs with MAF>0,05:  3242
	- Count and frequency of 0 (0/0), 1 (0/1) and 2 (1/1) for all variants:
		     0_count_obs  1_count_obs  2_count_obs
	rs555599246          150            0            0
	rs534396879          150            0            0
	rs116742944          138           12            0
	rs138086069          149            1            0
	rs535042364          150            0            0
	- SNPs with MAF>0.5, observed and expected counts:
		  rsID  alt_allele_frequency  ...  1_count_obs  2_count_obs
	6    rs6139892              0.836667  ...           29          111
	9    rs6053830              0.843333  ...           27          113
	12  rs62205718              0.073333  ...           22            0
	29   rs6085351              0.113333  ...           30            2
	38   rs6038330              0.113333  ...           30            2

	[5 rows x 9 columns]


	=======================================
	Q3. LD
	- Calculate pairwise D, D' and R2 of the first 100 SNPs
	.......... 500 pairs of SNPs processed
	.......... 1000 pairs of SNPs processed
	.......... 1500 pairs of SNPs processed
	.......... 2000 pairs of SNPs processed
	.......... 2500 pairs of SNPs processed
	.......... 3000 pairs of SNPs processed
	.......... 3500 pairs of SNPs processed
	.......... 4000 pairs of SNPs processed
	.......... 4500 pairs of SNPs processed
	.......... 5000 pairs of SNPs processed
	.
	In total: 5050 pairs of SNPs processed


	=======================================
	Q4. PCA
	- Replace 0/0, 0/1, 1/1 with 0, 1, 2 in a new DataFrame (df_additive_genotype):
		ID status  rs555599246  rs534396879  ...  rs577761377  rs6140157  rs185180844  rs527420293
	0  NA06984   case            0            0  ...            0          2            0            0
	1  NA06986   case            0            0  ...            0          2            0            0
	2  NA06989   case            0            0  ...            0          2            0            0
	3  NA07000   case            0            0  ...            0          2            0            0
	4  NA07037   case            0            0  ...            0          2            0            0

	[5 rows x 28746 columns]
	- Perform PCA, % variance explained by the first 10 PCs are:
		16.81%	7.31%	4.98%	3.62%	3.01%	2.63%	2.28%	2.21%	2.09%	1.95%

	# Run finished in 4.91 minutes
	```
