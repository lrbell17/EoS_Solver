# LB_EoS_Solver

This repository contains 2 python programs tovsolver_k2_v2.py and unit_conversion.py, as well as several Equations of States (EoS) for neutron stars. 

**tovsolver_k2_v2.py --** solves the Tolmann-Oppenheimer-Volkoff (TOV) Equation using the RK4 Method for a spectrum of stars for a given EoS. 

**Inputs:**
  - **EoS file (.dat):**
    - 3 Columns: Pressure (MeV/fm^3) | Central Energy Density (MeV/fm^3) | Baryon Density (1/fm^3)
  - **Initial Central Energy Density (e0) and Number of Stars (N):**
    - solves TOV Equation for e0 += 10 N times
    - For the program to behave properly e0 must be in the right range of values (typically between 100 and 2500 MeV/fm^3)
    - Examples:  
      - APR: e0 = 150, N = 150
      - SLY: e0 = 200, N = 110

**Outputs:** the user can plot any two values from the options listed below. This is an interactive plot that allows the user to click any point on the curve to output all computed values as well as the Mass/Radius curve for the single star at that position. 
  - Mass
  - Radius
  - Apsidal Constant (k2 or the 2nd Love Number)
  - Compactness 
  - Pressure
  - Central Energy Density
  - Tidal Polarizability (lambda)


**unit_conversion.py --** of the EoS files provided above, only APR (eosaapr.dat) and SLY (eossly.dat) are in the correct format. 

The rest have 4 Columns: Central Energy Density (g/cc) | Pressure (dynes/cm^2) | Enthalpy | Baryon Density (1/cm^3), which are not compatible with **tovsolver_k2_v2.py**. 

Thus, **unit_conversion.py** takes in an EoS (as described above) as an input and outputs an array with the correct format and units.


