# Influence of Databases

Files to reproduce results of results on https://arxiv.org/abs/2104.06099

For the Neural Network implementation please check: [PhysNet](https://github.com/MMunibas/PhysNet)

For the databases used: 
- [QM9](https://figshare.com/collections/Quantum_chemistry_structures_and_properties_of_134_kilo_molecules/978904)
- [PC9](https://figshare.com/articles/dataset/PC9_data_zip/9033977)
- [ANI-1E](https://zenodo.org/record/4680953)
- [ANI-1](https://figshare.com/collections/_/3846712)
- [ANI-1x](https://doi.org/10.6084/m9.figshare.c.4712477)
- [Tautobase-Smiles](https://github.com/WahlOya/Tautobase)
- [Tautobase-QC calculations](https://zenodo.org/record/4680972)

## Organization of the Files

- Databases: Contain the results for each database on .csv files. The structure of the file contains an id, Smiles for Isomer A, Smiles for Isomer B, number of heavy atoms(natoms), Energy of Isomer A from DFT calculation (Ea_DFT), Energy of Isomer A from NN calculation (Ea_NN), Difference for Isomer A (Diff_MolA), Energy of Isomer B from DFT calculation (Eb_DFT), Energy of Isomer B from NN calculation (Eb_NN), Difference for Isomer B (Diff_MolB),Tautomerization energy at DFT level(Diff_abi), Tautomerization energy with NN(Diff_NN), Error on Tautomerization Energy(Diff_NN_abi)
-  Kullback-Leibler (KL) divergence analysis: See inside the folder for more details.
- Scripts for plot: Scripts to reproduce the plots on the manuscript
- T-MAP: Eric
- Amons Analysis: Eric

