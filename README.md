# A trio of genes in germinal center B cells mediates the pathogenesis of blood cancers

Wei-Quan Fang\*, Walt Yu-Le Wu†, Ming-Jing Hwang\*
 
Institute of Biomedical Sciences, Academia Sinica, Taipei, Taiwan, R.O.C.

\*Corresponding authors. 

†Present address: European Molecular Biology Laboratory, Heidelberg, Germany.

## Introduction of the repository

This repository contains the raw data and matlab codes to reproduce the results in the manusript entitled: 
```python
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

```javascript
[Step1] Download the two files, CODEinGAN and DATAinGAN, and put them in C:\

[Step2] Open matlab program window and type *cd 'C:\CODEinGAN'* in it

[Step3] Type *GANreproduce* in matlab window and select dialog boxes for reproducing
```

Please see '''GANreproduce Quick Start Guide.pdf''' for more details.

## Interpretation of data files DATAinGAN

`ltaDS10gnD1`

`ltaDS10gnD2`

`ltaDS10gnD3`

`ltaDS10gnD4`

`ltaDS10gnD60`

`ltaUS5gnD1`

`ltaUS5gnD2`

`ltaUS5gnD3`

`ltaUS5gnD4`

`ltaUS5gnD60`


## Interpretation of code files CODEinGAN

`logrankLF.m`

`LT3X2Ana.m`

## Citation and Contact

A manuscript detailing our work has been submitted and the citation will be provided later on. Any question or suggestion please feel free to contact [WQ Fang](mailto:deleapoli@gmail.com) or [MJ Hwang](mailto:mjhwang@ibms.sinica.edu.tw).
