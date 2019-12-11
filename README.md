# eSTARR code base
eSTARR-seq and HiDRA Analysis scripts accompanying Tippens & Liang et al (https://www.biorxiv.org/content/10.1101/818849v1.full)


The get_datafiles.sh will need to be run first to obtain larger datafiles from public repositories.

Also, strand-specific HiDRA read counts (hg19) must be downloaded [here](https://drive.google.com/open?id=1fdxk-D_2M-TVV3C-IXw5oNAsZeCenRU_).

Key analyses are found in the following files:

* HiDRA_voom.ipynb: Exemplary read count analysis of strand-specific HiDRA data using voom.
* HiDRA_annot.ipynb: Exemplary analyses of HiDRA fragments around GRO-cap TSSs.
* eSTARR_EnhCalls.ipynb: Exemplary read count analysis of eSTARR data using voom. Estimate strand and size bias.
* Enh_Fusions_Coop.ipynb: Exemplary analysis of enhancer clusters from eSTARR data.
