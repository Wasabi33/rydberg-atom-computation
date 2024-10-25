rydberg-atom-computation
(2023 summer UROP project)
==========================


This is a computational project about rydberg atoms. A python package called ARC has been used to calculate different properties and parameters of rydberg atoms. We invesigated how to prepare rydberg states using both one-photon and two-photon processes and found the optimal parameters e.g. laser detuning to achieve the best transfer efficiency. We developed a more efficent computational method to solve the linear optical bloch equations and compared with the more traditional way of solving Born-Markov master equation directly. Next, we investigated numerically the effect of adiabatic rapid passage(ARP) on transferring the atoms to the circular rydberg states (with m=l): the ARP map has been calculated. Due to the break of symmetry by applying a z-direction electric field, l and n are no longer good quantum numbers. We therefore switched to a new coordinate system to use parabolic quantum numbers. In the new coordinate system, the stark map was calculated and a 'diamond manifold' was seen, which is the key to transfer the population to high m number circular rydberg states through the ARP process. Finally we calculated and plotted the wavefunction of rydberg atoms during the process of ARP and could visualise the specific symmetry in the parabolic coordinate system.

More details can be found in this presentation `Preparation_of_high_principal_quantum_number_Rydberg_states.pdf`

References are in the reference folder.

----------
Environment
----------
The project uses the (Alkaline Rydberg Calculator) ARC package, which is a python package for calculating properties of alkali-metal Rydberg atoms. More details can be found on: https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/

-------
Authors
-------

[Chen Lu](https://github.com/Wasabi33), Supervised by: [Prof Mike Tarbutt](https://www.imperial.ac.uk/people/m.tarbutt)

This is an Undergraduate Research Opportunities Programme (UROP) project at Imperial College London, 2023 summer.