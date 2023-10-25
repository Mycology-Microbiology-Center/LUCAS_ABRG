[![DOI](https://zenodo.org/badge/708775478.svg)](https://zenodo.org/doi/10.5281/zenodo.10039980)

# A trait-based approach unveils mechanisms shaping antibiotic-related microbial genetic machinery across European soils <img src='images/MMC_logo.png' align="right" height="100" />

This repository houses the source code related to the study of antibiotic synthesis and resistance genes, obtained using high-throughput sequencing to analyze environmental DNA from 658 metagenomic samples of topsoil, collected from various regions across Europe.  

## Materials

The foundation for the study is based on the European Commission's [LUCAS Soil survey](https://esdac.jrc.ec.europa.eu/projects/lucas). 
These samples were collected across the European Union using a standardized protocol.  

Distribution of sampling sites across environments:  
<p align="middle">
  <img src="images/sampling_map.jpg" width="600" title="Dataset"/>
</p>


## Repository structure

- **Main directory**: contains the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow files and R-scripts  
- **`scripts` subdirectory**: contains auxiliary scripts utilized in the workflow  
- **`envs` subdirectory**: holds the specifications for [conda](https://docs.conda.io/en/latest/) environments  
- **`Singularity_Containers` subdirectory**: stores the [Singularity](https://sylabs.io/singularity/) container definition files  


## Citation

Dulya O, Mikryukov V, Shchepkin DV, Pent M, Tamm H, Guazzini M, Panagos P, Jones A, Orgiazzi A, Marroni F, Bahram M, Tedersoo L. (2023) 
A trait-based approach unveils mechanisms shaping antibiotic-related microbial genetic machinery across European soils // [DOI:10.5281/zenodo.10039980](https://zenodo.org/doi/10.5281/zenodo.10039980)


