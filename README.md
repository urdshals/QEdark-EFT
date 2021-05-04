# QEdark-EFT

QEdark-EFT is an extention of the [QEdark code](https://github.com/adrian-soto/QEdark_repo) previously used by [Essig et. al.](https://arxiv.org/abs/1509.01598) to [calculate the crystal form factor](http://ddldm.physics.sunysb.edu/ddlDM/). Both QEdark-EFT and QEdark run on [QuantumESPRESSO](https://www.quantum-espresso.org), and requires [QuantumESPRESSO v.5.1.2](https://github.com/QEF/q-e/releases/tag/qe-5.1.2) to run [1](https://iopscience.iop.org/article/10.1088/0953-8984/21/39/395502) [2](https://aip.scitation.org/doi/10.1063/5.0005082).

As described in this publication, QEdark-EFT extends the calculation of the crystal form factor to calculate crystal responses required to calculate the excitation rate in Silicon and Germanium caused by general nonrelativistic dark matter.

## Summary

- In the paper we show that any nonrelativistic dark matter electron interaction in Silicon and Germanium can be expressed in terms of 5 products of dark matter and crystal responses. 
- These crystal responses are expressed in terms of electron wave-function overlap integrals, which are evaluated by QEdark-EFT using the DFT tools from QuantumESPRESSO.
- The code written in FORTRAN90 and uses [openMP](https://www.openmp.org) for parallelisation.

<details><summary>Repository content</summary>
<p>

The included folders are:

- *run/*: Contains the input files used for the final runs which generated the crystal responses shown in the paper.
- *test/*: Contains the input files for test runs intended to finish in a few minutes on a laptop.
- *QEdark-EFT/*: Contains the files that need to be copied into the QuantumESPRESSO folder. 

</p>
</details>

## Installing and running QEdark-EFT


<details><summary>1. Downlad & Installation</summary>
<p>

1. Download and unzip [QuantumESPRESSO v.5.1.2](https://github.com/QEF/q-e/releases/tag/qe-5.1.2)
2. Copy the files in in *QEdark-EFT/* into the QuantumESPRESSO directory. Some of the files will replace files that are already there with the same name.
3. While in the QuantumESPRESSO main directory execute ``` make clean ```, then ``` ./configure --disable-parallel --enable-openmp ```, then ``` make pw ```. The first two commands are just needed the first time the code is installed, the latter is needed every time you have made changes to one of the files in the QuantumESPRESSO directory.
</p>
</details>
<details><summary>2. Usage</summary>
<p>

1. Copy *PW/src/pw.x* to the directory that contains the input files, for instance *test/Si/*
2. Set the number of threads and memory per thread. This is done by executing the commands  ``` export OMP_NUM_THREADS=n ``` and ``` export OMP_STACKSIZE=mM ```, where n is the number of threads and m is the memory in MB used per thread. For the test run one can use ``` export OMP_NUM_THREADS=8 ``` and ``` export OMP_STACKSIZE=500M ```.
3. Run the code from the directory containing the input files by executing i.e. ``` ./pw.x <si.test.in> si.out ```. 
</p>
</details>

## Version History

- 04.05.2021: Release of version 0.1.0

## Authorship and Lisnece

<details><summary>Citation</summary>
<p>

If you decide to use this code, please cite the latest archived version and the paper.

</p>
</details>

<details><summary>Author & Contact</summary>
<p>

The author of QEdark-EFT is Einar Urdshals.

For questions, bug reports or other suggestions please contact Einar Urdshals (urdshals@chalmers.se).
</p>
</details>

<details><summary>License</summary>
<p>

This project is licensed under the MIT License - see the LICENSE file.

</p>
</details>
