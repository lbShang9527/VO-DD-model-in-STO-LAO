# VO-DD-model-in-STO-LAO
The Drift-Diffusion Modeling of Oxygen Vacancy Migration in STO/LAO. 

The N_d_0 file is the initial charge distribution file of VOs, which is a balanced distribution under zero voltage.
The dos file describes the relationship between the bottom electron concentration of the STO conduction band and the Fermi level.
The main.py is the main program. 

Four steps of input need to be adjusted during running. 
   1.  Set the bias voltage and position of the tip.
   2.  Read the initial concentration distribution file of VO.
   3.  Run the main iteration equation and set the number of steps.
   4.  Chose which output file to be printed.

There are three main output files.
The N_d file is the charge distribution file of VOs, whose unit is angstrom^-2.//
The N_e file is the charge distribution file of electrons, whose unit is angstrom^-2.//
The N_list file contains four columns of data, recording the total number of electrons in LAO, the total number of electrons in STO, the total number of VOs in LAO, and the total number of VOs in STO in each iteration.
