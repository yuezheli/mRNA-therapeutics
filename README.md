# How to use this repo

See readme files in each folder to learn more about the content. 

## General setup

### package management

1. Install pkgr software following the instructions [here](https://github.com/metrumresearchgroup/pkgr). 
2. Open the R project `mRNA.Rproj`. This allows you to work from within a self-contained project environment.
3. Install packages by typing in terminal: `pkgr install`. This command will look for the file `pkgr.yml` and install the packages listed. The specific package versions are imported from https://mpn.metworx.com/docs/.

### mrgsolve installation
For detailed instructions on mrgsolve installation and important dependencies, follow this link https://github.com/metrumresearchgroup/mrgsolve

Switch directory to folder-of-interest to reproduce all the figures. 

# Content of the folder

- README.md (this readme file)

- mRNA.Rproj (R project file)

- ```.gitignore``` (This file contains types of files that are not included in this GitHub repo)

- ```pkgr.yml``` (The project file lists the package version and source)

- Apgar2018: Implementation of model from [Apgar et al., 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391595/)

- Banks2003: Implementation of models from [Banks et al., 2003](https://www.nature.com/articles/3302076)

- Kagan2013: Implementation of models from [Kagan et al., 2014](https://pubmed.ncbi.nlm.nih.gov/23793994/). Note this models a liposomal distribution of small molecule drug, not really a gene therapy. 

- Mihaila2017: Implementation of model from [Mihaila et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415968/)

- Varga2005: Implementation of model from both [Varga et al., 2001](https://pubmed.ncbi.nlm.nih.gov/11708880/) and [Varga et al., 2005](https://www.nature.com/articles/3302495)

# Gene therapy model comments

[Ledley and Ledley, 1994](https://pubmed.ncbi.nlm.nih.gov/7948130/) is one of the earliest model for gene therapy. This model focuses on naked DNA plasmid. The model includes DNA uptake from extracellular environment to cytosol, followed by transcription and translation, and protein secretion. The biggest caveat is that the model is not validated by experimental data. 

[Banks et al., 2003](https://www.nature.com/articles/3302076) focuses on plasmid DNA. It modeled plasmid DNA entering cytosol and into nucleus. The model was validated using fluorescent-labeled plasmids. Plasmid amount in cytosol and in nucleus were quantified based on fluorescence data. However, this model does not include any process related to protein production. The process for DNA uptake was also coarse. 

[Varga et al., 2005](https://www.nature.com/articles/3302495) proposed a more detailed mechanistic model for DNA uptake. This in vitro model included plasmid uptake, endosomal escape, nuclear import, plasmid unpacking in both cytosol and nucleus, and translation (though simplified). Parameters were obtaine from 3 sources: 

+ Fluorescence-based assays for cellular binding and uptake. 
+ Isolated DNA from cytosol and nuclear
+ Literature values

This model was fitted for multiple vectors/ transfection methods, including adenovirus 5.

[Mihaila et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415968/) modeled siRNA-mediated therapy in vitro. The model could be divided into 3 critical steps:

+ crossing plasma membrane, endosomal escape/unpackaging
+ loading of siRNA on to RISC (RNA induced silencing complex)
+ mRNA knockdown

Parameters used in this model was based on 4 different sources: 

+ LNP-siRNA complex uptake through cell-associated confocal microscopic imaging of labeled-siRNA
+ mRNA expression over time
+ binding assays for siRNA-RISC complex kinetics
+ literature values

The model identified the rate-limiting steps to be endosomal scaping and RISC complex formation.

[Apgar et al., 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391595/) modeled mRNA-LNP therapy. This model was used for first-in-human dose projections. The model contained 3 different parts: 

+ known mechanisms of LNP stability, biodistribution
+ attachment to hepatocytes, endocytic uptake, intracellular unpacking of LNPs
+ protein expression of UGT1A1, and enzymatic modification of bilirubin leading to its increased clearance.

The data used to validate the model was obtained from Gunn rats, including

+ Plasma levels of intact LNP
+ Liver mRNA levels
+ Total bilirubin levels in Plasma
+ Protein Expression in Liver

The model was humanized by

+ Bilirubin synthesis and degradation rates are updated
+ Biophysical interactions or cellular processes are expected to be preserved between humans and rats
+ Protein degradation rate and rate of bilirubin catalysis were updated
+ Volumes were updated to reflect human physiology
