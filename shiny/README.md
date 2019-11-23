### Overview
This tool allows users to estimate the power of their sequencing experiments for CNV detection.  
The formulae used for power calculations are described in: "https://www.biorxiv.org/content/10.1101/413690v1"  
If you find this tool useful, please consider cite our paper.  

### Parameters specified by user
F: the portion of the sample that contains the CNA or CNV  
alpha: the significance level  
L: the length of the CNA or CNV  
N: the ploidy of the CNA or CNV (normal diploid regions have N=2)  
W: the size of the window  
l: average sequencing read length  
D: haploid sequencing depth (for a diploid dataset of 30X total coverage, the haploid coverage is 15)  
theta: the variance inflation factor, equal to phi*N*D*W/l  
phi: the variance inflation factor

### Contact
This website is maintained by Jun Li's Lab.  
Please contact hanyou@umich.edu and sramdas@umich.edu for questions.

### Source code
GitHub: [github.com/shwetaramdas/cnvpoweranalysis](http://github.com/shwetaramdas/cnvpoweranalysis)  
Web App: [cnvpowercalculator.junzli.com](http://cnvpowercalculator.junzli.com)  
