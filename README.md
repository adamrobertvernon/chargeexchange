# Simulation of atomic populations of ions following charge exchange (1-100 keV)
Python and mathematica scripts for simulation of atom relative populations following charge-exchange

If you do not know what this program is for please check these papers:
Vernon, A.R. et al. “Simulation of the Relative Atomic Populations of Elements 1 ≤ Z ≤89 Following Charge Exchange Tested with Collinear Resonance Ionization Spectroscopy of Indium.” Spectrochimica Acta Part B: Atomic Spectroscopy 153 (March): 61–83. https://doi.org/10.1016/J.SAB.2019.02.001.

Ryder, C.a., K. Minamisono, H.B. Asberry, B. Isherwood, P.F. Mantica, a. Miller, D.M. Rossi, and R. Strum. 2015. “Population Distribution Subsequent to Charge Exchange of 29.85keV Ni+ on Sodium Vapor.” Spectrochimica Acta Part B: Atomic Spectroscopy 113: 16–21. https://doi.org/10.1016/j.sab.2015.08.004.

Rapp, Donald, and W. E. Francis. 1962. “Charge Exchange between Gaseous Ions and Atoms.” The Journal of Chemical Physics 37 (11): 2631. https://doi.org/10.1063/1.1733066.

# Setup
Requires Mathematica >10 
Tested only on Linux so far

## Create a script to run Mathematica commands from Python (http://sapiensgarou.blogspot.com/2012/06/how-to-run-mathematica-functions-on.html):

Create a script named runMath with the content:

    #!/usr/local/bin/MathematicaScript -script

    value=ToExpression[$ScriptCommandLine[[2]]];

    Print[value];

Give execution privilege to the file:

    sudo chmod +x runMath

Move the file to the execution path:

    sudo mv runMath /usr/bin/
    
If you change the location of the script you must also update "RF_mathematica.py"

## Running the simulations script 
### Modify the parameters in sim_list_all.csv file for the atomic system/s you are interested in.

+ Transitions - change this to e.g. [0.0, 40636.98] for an lower and upper state, if you want the satellite peak contributions to also be calculated, else put it to False
+ collinear? - True for collinear geometry, False for anti-collinear. This is used for the satellite peak calculations.
+ T_B - beam energy in eV
+ dist - distance of atom flight for the simulated final atomic population
+ time_steps - number of step in time for the population distribution, the actual time periods will depend on the beam velocity. More is better. The number for accurate populations and convergence will depend on the transitions of the atomic system.
+ skip - Change skip to 0 to let that row be calculated, you can include many rows for systematic calculations.
+ database - nist by default. You can also use kurucz, see the example .kurucz file  format in the '/datatables' folder, you will need to download the data from kurucz for your atomic system
+ state - leave as 1 for neutralisation of 1+ ions to their neutral state. It is not recommended to use these calculations for 2+ to 1+ but this can be performed by downloading the appropriate element data. See the "FII..." files for examples and change 'state' to 2

### run 'run_popdist_sims.py'
The results will be saved in the '/results' folder.
Modify the program as you need.
