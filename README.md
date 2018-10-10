# A trio of genes in germinal center B cells mediates the pathogenesis of blood cancers

Wei-Quan Fang\*, Walt Yu-Le Wu?, Ming-Jing Hwang\*
 
Institute of Biomedical Sciences, Academia Sinica, Taipei, Taiwan, R.O.C.

\*Corresponding authors. 

?Present address: European Molecular Biology Laboratory, Heidelberg, Germany.

## Introduction of the repository

This repository contains the raw data and matlab codes to reproduce the results in the manusript entitled: 
```
A trio of genes in germinal center B cells mediates the pathogenesis of blood cancers.
```

In this work, we analyzed lymphoma-related gene expression and clinical data and identified a regulatory motif of germinal center B cell genes. This motif involves a trio of genes and the regulation of two oncogenes and plays crucial roles in blood cancers. There are three significant advances:

```
1. This motif was largely consistent across five different transcriptomic datasets.

2. The inference from this motif provided a possible resolution to a clinically observed controversy on the effect of BACH2 on lymphoma.

3. A common cause and a possible mechanism involving retro-differentiation for pan blood cancer was suggested.
```

According to these findings, we can conduct multiple blood cancer prognosis and uncover tumor supressor functions of BACH2 which are much better than previous thought.

## Reproducing step

Note that an installed matlab program is a must for the following steps.

* **_Step1_** Download the two files, `CODEinGAN` and `DATAinGAN`, and put them in C:\

* **_Step2_** Open matlab program window and type `cd 'C:\CODEinGAN'` in it

* **_Step3_** Type `GANreproduce` in matlab window and select dialog boxes for reproducing

Please see  `GANreproduce Quick Start Guide.pdf`  for more details of the GAN (gene association network) analysis.

## Interpretation of data files DATAinGAN

`realDS10gn` is about down-stream genes of GAN used in GSE60 for real data analysis.

`realUS5gn`  is about up-stream genes of GAN used in GSE60 for real data analysis.

`ltaDS10gnD1` is about down-stream genes of GAN used in GSE4475 for linear trend analysis.

`ltaDS10gnD2` is about down-stream genes of GAN used in GSE32918 for linear trend analysis.

`ltaDS10gnD3` is about down-stream genes of GAN used in GSE31312 for linear trend analysis.

`ltaDS10gnD4` is about down-stream genes of GAN used in GSE10846 for linear trend analysis.

`ltaDS10gnD60` is about down-stream genes of GAN used in GSE60 for linear trend analysis.

`ltaUS5gnD1` is about up-stream genes of GAN used in GSE4475 for linear trend analysis.

`ltaUS5gnD2` is about up-stream genes of GAN used in GSE32918 for linear trend analysis.

`ltaUS5gnD3` is about up-stream genes of GAN used in GSE31312 for linear trend analysis.

`ltaUS5gnD4` is about up-stream genes of GAN used in GSE10846 for linear trend analysis.

`ltaUS5gnD60` is about up-stream genes of GAN used in GSE60 for linear trend analysis.

`sDLBCL1` contains clinical survival and trio's expression data in GSE4475 (DLBLC).

`sDLBCL2` contains clinical survival and trio's expression data in GSE32918 (DLBLC).

`sDLBCL3` contains clinical survival and trio's expression data in GSE31312 (DLBLC).

`sDLBCL4` contains clinical survival and trio's expression data in GSE10846 (DLBLC).

`sOSBC2a` contains clinical survival and trio's expression data in GSE16131 (FL).

`sOSBC2b` contains clinical survival and trio's expression data in GSE66166 (FL).

`sOSBC3a` contains clinical survival and trio's expression data in GSE9782 (MM).

`sOSBC3b` contains clinical survival and trio's expression data in GSE24080 (MM).

`sOSBC4` contains clinical survival and trio's expression data in GSE22762 (CLL).

`sOSBC5` contains clinical survival and trio's expression data in Rosenwald (MCL).

`sOSBC6a` contains clinical survival and trio's expression data in TCGA (AML).

`sOSBC6b` contains clinical survival and trio's expression data in GSE12417 (AML).

`sOSBC7` contains clinical survival and trio's expression data in TARGET (ALL).


## Interpretation of code files CODEinGAN

`GANreproduce` reproduces the results in the manuscript via **_User Friendly Dialog Boxes_**.

`logrankLF.m` generates log-rank test by package of [Fan Lin](https://www.mathworks.com/matlabcentral/fileexchange/20388).

`LT3X2Ana.m` generates consistent linear trend analysis for the motif. 

`realAna` generates real data analysis for the AWTE-induced GAN.

`regAWTE` returns regression model fitting by AWTE method.

`regLasso1` returns regression model fitting by LASSO method.

`rrHKB` returns regression model fitting by RR method.

`simuAna` generates simulated data analysis for four different methods.

`simuAwAna` generates simulated data analysis for robustness of AWTE method.

`survBAna1` generates clinical controversy for BACH2's role.

`survBAna2` generates BACH2 higher and survival better in patients with lower SPIB expression.

`survBAna3` generates BACH2 higher and survival worse in patients with higher SPIB expression.

`survTAna1` produces distinguishable survival prognosis for the trio in DLBCL.

`survTAna2` produces distinguishable survival prognosis for the trio in FL.

`survTAna3` produces distinguishable survival prognosis for the trio in MM.

`survTAna4` produces distinguishable survival prognosis for the trio in CLL.

`survTAna5` produces distinguishable survival prognosis for the trio in MCL.

`survTAna6` produces distinguishable survival prognosis for the trio in AML.

`survTAna7` produces distinguishable survival prognosis for the trio in ALL.


## Citation and Contact

A manuscript detailing our work has been submitted and the citation will be provided later on. Any question or suggestion please feel free to contact [WQ Fang](mailto:deleapoli@gmail.com) or [MJ Hwang](mailto:mjhwang@ibms.sinica.edu.tw).
