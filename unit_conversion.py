'''
Input: .dat file with columns Central Energy Density (g/cc),
Pressure (dynes/cm^2), Enthalpy (which can be ignored), and Baryon Density
(1/cm^3).

Outputs: columns Pressure (MeV/fm^3), Central Energy Density (MeV/fm^3),
and Baryon Density(1/fm^3), which is compatible with tovsolver_k2_v2.py
'''



import numpy as np

L = []
while (len(L) ==0 or len(L[0,:]) !=4):
    try:
        EoS_file = input('\nPlease give the name of the EoS file (eosaapr.dat):\t')
        while(EoS_file[len(EoS_file)-4:len(EoS_file)] != '.dat'):
            print('\n\tERROR -- Please enter .dat file\n')
            EoS_file = input('Please give the name of the EoS file (eosaapr.dat):\t')
    
        L = np.loadtxt(EoS_file)   #--Reading the EoS.
        if len(L[0,:]) !=4:
            print('\n\tERROR -- Incorrect dimensions\n')
        
    except Exception:
        print('\n\tERROR -- File name is not correct!\n')
        pass
    

M0 = []
M1 = []
M2 = []

k = 1.6022*10**33       #-- 1 Mev/fm^3 = k*[dynes/cm^2]
l = 1.7827*10**12       #-- 1 MeV/fm^3 = l*[g/cc]
m = 10**39              #-- 1/fm^3 = m*[1/cm^3]

for i in range(0,len(L)-1):
    
    M0.append(L[i,1]/k)
    M1.append(L[i,0]/l)
    M2.append(L[i,3]/m)
    print(M0[i], M1[i], M2[i])   
