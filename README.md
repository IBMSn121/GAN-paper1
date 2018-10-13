# A trio of genes in germinal center B cells mediates the pathogenesis of blood cancers

Wei-Quan Fang\*, Walt Yu-Le Wu?, Ming-Jing Hwang\*
 
Institute of Biomedical Sciences, Academia Sinica, Taipei, Taiwan, R.O.C.

\*Corresponding authors. 

?Present address: European Molecular Biology Laboratory, Heidelberg, Germany.

## Introduction of the repository

This repository contains the raw data and matlab codes to reproduce the results reported in the manuscript: 
```
A trio of genes in germinal center B cells mediates the pathogenesis of blood cancers.
```

In this work, we analyzed lymphoma-related gene expression and clinical data and identified a regulatory motif of germinal center B cell genes. This motif involves a trio of genes and the regulation of two oncogenes. The motif was shown to play crucial roles in blood cancers. 

## Reproducing step

Note that an installed matlab program is required to use the GANreproduce Package in the following steps.

* **_Step1_** Download the two files, `DATAinGAN` and `CODEinGAN`, and put them in C:\

* **_Step2_** Open matlab program window and type `cd 'C:\CODEinGAN'` in it

* **_Step3_** Type `GANreproduce` in the matlab window and select dialog boxes to reproduce the results

Please see `GANreproduce Quick Start Guide.pdf` for more details and examples of reproducing the GAN (gene association network) related analyses.

## Interpretation of data file DATAinGAN

`realDS10gn.txt` is about down-stream genes of GAN used in GSE60 for analysis of the actual lymphoma data.

`realUS5gn.txt`  is about up-stream genes of GAN used in GSE60 for analysis of the actual lymphoma data.

`ltaDS10gnD1.txt` is about down-stream genes of GAN used in GSE4475 for linear trend analysis.

`ltaDS10gnD2.txt` is about down-stream genes of GAN used in GSE32918 for linear trend analysis.

`ltaDS10gnD3.txt` is about down-stream genes of GAN used in GSE31312 for linear trend analysis.

`ltaDS10gnD4.txt` is about down-stream genes of GAN used in GSE10846 for linear trend analysis.

`ltaDS10gnD60.txt` is about down-stream genes of GAN used in GSE60 for linear trend analysis.

`ltaUS5gnD1.txt` is about up-stream genes of GAN used in GSE4475 for linear trend analysis.

`ltaUS5gnD2.txt` is about up-stream genes of GAN used in GSE32918 for linear trend analysis.

`ltaUS5gnD3.txt` is about up-stream genes of GAN used in GSE31312 for linear trend analysis.

`ltaUS5gnD4.txt` is about up-stream genes of GAN used in GSE10846 for linear trend analysis.

`ltaUS5gnD60.txt` is about up-stream genes of GAN used in GSE60 for linear trend analysis.

`sDLBCL1.txt` contains clinical survival and trio's expression data in GSE4475 (DLBLC).

`sDLBCL2.txt` contains clinical survival and trio's expression data in GSE32918 (DLBLC).

`sDLBCL3.txt` contains clinical survival and trio's expression data in GSE31312 (DLBLC).

`sDLBCL4.txt` contains clinical survival and trio's expression data in GSE10846 (DLBLC).

`sOSBC2a.txt` contains clinical survival and trio's expression data in GSE16131 (FL).

`sOSBC2b.txt` contains clinical survival and trio's expression data in GSE66166 (FL).

`sOSBC3a.txt` contains clinical survival and trio's expression data in GSE9782 (MM).

`sOSBC3b.txt` contains clinical survival and trio's expression data in GSE24080 (MM).

`sOSBC4.txt` contains clinical survival and trio's expression data in GSE22762 (CLL).

`sOSBC5.txt` contains clinical survival and trio's expression data in Rosenwald (MCL).

`sOSBC6a.txt` contains clinical survival and trio's expression data in TCGA (AML).

`sOSBC6b.txt` contains clinical survival and trio's expression data in GSE12417 (AML).

`sOSBC7.txt` contains clinical survival and trio's expression data in TARGET (ALL).


## Interpretation of code file CODEinGAN

`GANreproduce.m` reproduces the results in the manuscript via **_user friendly dialog boxes_**.

`logrankLF.m` generates log-rank test by package of [Fan Lin](https://www.mathworks.com/matlabcentral/fileexchange/20388).

`LT3X2Ana.m` generates consistent linear trend analysis for the motif. 

`realAna.m` generates analysis of the actual lymphoma data for the AWTE-induced GAN.

`regAWTE.m` returns regression model fitting by AWTE method.

`regLasso1.m` returns regression model fitting by LASSO method.

`rrHKB.m` returns regression model fitting by RR method.

`simuAna.m` generates simulated data analysis for four different methods.

`simuAwAna.m` generates simulated data analysis for robustness of AWTE method.

`survBAna1.m` generates clinical controversy for BACH2's role.

`survBAna2.m` generates results of BACH2 higher and survival better in patients with lower SPIB expression.

`survBAna3.m` generates results of BACH2 higher and survival worse in patients with higher SPIB expression.

`survTAna1.m` produces distinguishable survival prognosis for the trio in DLBCL.

`survTAna2.m` produces distinguishable survival prognosis for the trio in FL.

`survTAna3.m` produces distinguishable survival prognosis for the trio in MM.

`survTAna4.m` produces distinguishable survival prognosis for the trio in CLL.

`survTAna5.m` produces distinguishable survival prognosis for the trio in MCL.

`survTAna6.m` produces distinguishable survival prognosis for the trio in AML.

`survTAna7.m` produces distinguishable survival prognosis for the trio in ALL.


## Citation and Contact

A manuscript detailing our work has been submitted and the citation will be provided later. Any question or suggestion please feel free to contact [WQ Fang](mailto:deleapoli@gmail.com) or [MJ Hwang](mailto:mjhwang@ibms.sinica.edu.tw).
