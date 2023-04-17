[![License: GPL 
v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![arXiv](https://img.shields.io/badge/arXiv-2303.15497-B31B1B.svg)](https://arxiv.org/abs/2303.15497)
[![arXiv](https://img.shields.io/badge/arXiv-2303.15509-B31B1B.svg)](https://arxiv.org/abs/2303.15509)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7836577.svg)](https://doi.org/10.5281/zenodo.7836577)


# QEdark-EFT for graphene with Monte Carlo Integrator

QEdark-EFT for graphene is an extention of the [QEdark code](https://github.com/adrian-soto/QEdark_repo) previously used by [Essig et. al.](https://arxiv.org/abs/1509.01598) to [calculate the crystal form factor](http://ddldm.physics.sunysb.edu/ddlDM/). QEdark-EFT for graphene runs on [QuantumESPRESSO](https://www.quantum-espresso.org), and requires [QuantumESPRESSO v.6.4.1](https://github.com/QEF/q-e/releases/tag/qe-6.4.1) to run [[1]](https://iopscience.iop.org/article/10.1088/0953-8984/21/39/395502) [[2]](https://aip.scitation.org/doi/10.1063/5.0005082).

As described in the publications, QEdark-EFT for graphene calculates the single material response function needed to calculate the rate of electron ejections caused by general nonrelativistic dark matter. The Monte Carlo Integrator uses this material response function to obtain observable ejection rates.

## Summary

- In the papers we show that any nonrelativistic dark matter induced electron ionisation in graphene can be expressed in terms of a single product of a dark matter and a material response function. 
- The material response function is expressed in terms of the momentum space electron density, which is evaluated by QEdark-EFT for graphene using the DFT tools from QuantumESPRESSO.
- QEdark-EFT for graphene is written in FORTRAN90 and uses [openMP](https://www.openmp.org) for parallelisation. For running the code IFX and openMP has been used. The Monte Carlo Integrator is written in Python and uses Joblib, Math, Numpy, Os, Scipy, Sys and Time.

<details><summary>Repository content</summary>
<p>

The included folders are:

- *Monte Carlo Integrator/*: Contains the QEdark-EFT input files required to compute the material response function, as well as the Python files that load the material response function and computes the observable ejection rates.
- *QEdark-EFT/*: Contains the files that need to be copied into the QuantumESPRESSO folder. 

</p>
</details>

## Installing and running QEdark-EFT


<details><summary>1. Downlad & Installation</summary>
<p>

1. Download and unzip [QuantumESPRESSO v.6.4.1](https://github.com/QEF/q-e/releases/tag/qe-6.4.1)
2. Copy the files in *QEdark-EFT/* into the QuantumESPRESSO directory. Some of the files will replace files that are already there with the same name.
3. While in the QuantumESPRESSO main directory execute ``` make clean ```, then ``` ./configure --enable-openmp ```, then ``` make pw ```. The first two commands are just needed the first time the code is installed, the latter is needed every time you have made changes to one of the files in the QuantumESPRESSO directory. 
</p>
</details>
<details><summary>2. Usage</summary>
<p>

1. Set up *Monte Carlo Integrator/data/high_ecut/run_scf.sh* and *Monte Carlo Integrator/data/low_ecut/run_scf.sh* to run on your cluster and submit them using sbatch. When they have run, there should be files named *E.dat*, *kG.dat* and *W.dat* in *Monte Carlo Integrator/data/high_ecut/* and *Monte Carlo Integrator/data/low_ecut/*.
2. Set up the file *Monte Carlo Integrator/run.sh* to run on your cluster. Select the DM masses you wish to compute the rates for in *Monte Carlo Integrator/all_rate_computation.py*, and set the interaction type and choose between CNTs and graphene in *Monte Carlo Integrator/run_all_rate.py*.
3. Run the Monte Carlo code by running the Python script *run_all_rate.py* in the *Monte Carlo Integrator/* directory, i.e. ``` python3 run_all_rate.py ```. 
</p>
</details>

## Version History

- 17.04.2023: Release of version 0.1.0

## Authorship and License

<details><summary>Citation</summary>
<p>

If you decide to use this code, please cite the latest archived version and the papers.

> Urdshals, E, & Matas, M. (2023). QEdark-EFT-Graphene (0.1.0). [[Zenodo.7836577]](https://doi.org/10.5281/zenodo.7836577)

> Catena, R., Emken, T., Matas, M, Spaldin, N.A., Urdshals, E. , 2023,  **Direct searches for general dark matter-electron interactions with graphene detectors: Part I. Electronic structure calculations**, [[arXiv:2303.15497]](https://arxiv.org/abs/2303.15497).
  
> Catena, R., Emken, T., Matas, M, Spaldin, N.A., Urdshals, E. , 2023,  **Direct searches for general dark matter-electron interactions with graphene detectors: Part II. Sensitivity studies**, [[arXiv:2303.15509]](https://arxiv.org/abs/2303.15509).

</p>
</details>

<details><summary>Author & Contact</summary>
<p>

The author of QEdark-EFT is Einar Urdshals and Marek Matas.

For questions, bug reports or other suggestions please contact Einar Urdshals (urdshals@chalmers.se).
</p>
</details>

<details><summary>License</summary>
<p>

This project is licensed under the [GPL License](http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt) - see the LICENSE file.

</p>
</details>
