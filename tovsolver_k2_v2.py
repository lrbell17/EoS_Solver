########################################################
###            Evaluating MR_curve & yR              ###
###                    SI Unit                       ###
########################################################

"""
This program can be used to plot the variation in Mass, Radius,
the apsidal constant k2, Pressure, Central Energy Density, Compactness,
and Tidal Polarizability for a given Equation of States (EoS). 

Runge-Kutta 4 is used to simultaneously solve for Mass, radius, and k2. 
The program can be run using the linux terminal. Python libraries:
'numpy', 'scipy.constants', 'scipy', 'matplotlib.pyplot',
'matplotlib.ticker', and 'mplcursors' are required to use the program.

After inputting initial values, wait 30-60s and follow instructions.
        
"""

import math
import numpy as np
from numpy import pi
import scipy.constants as const
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import mplcursors



M_sun = 1.98847*10**30
G = const.G
c = const.c
k = G/c**2                #--this value is in N/kg^2.s^2 ~ m/kg
fact1 = G/c**4
log = math.log   


#--Function for float range--#
def frange (start, stop, step):
    i = start
    while i<stop:
        yield i
        i += step


       
#--Plotting values for all stars--#
def plotting (X, Y, mass, radius, P, E, k2, C, Lam, RADIUS, MASS, marker, xname, yname, title):    
    rvals = []
    mvals = []
    
    fig, ax = plt.subplots(1)
    graph = ax.plot(X, Y, marker)
    ax.set_xlabel(xname, fontsize=15)
    ax.set_ylabel(yname, fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.grid(True, which='major')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=15, length=10)
    ax.tick_params(axis='both', which='minor', length=5)

    cursor = mplcursors.cursor(graph)   # Allows user to interact w/ graph and get values
    @cursor.connect("add")
    def return_values(sel):
        tup1 = cursor.selections   # Gets x&y values for selection
        xx, yy = tup1[0][1]
        
    #-- Getting index of closest values:
        idx = min(range(len(X)), key=lambda i: abs(X[i]-xx))  
        R_val = radius[idx]
        M_val = mass[idx]
        P_val = P[idx]
        E_val = E[idx]
        k2_val = k2[idx]
        C_val = C[idx]
        Lam_val = Lam[idx]
        rvals.append(RADIUS[idx])
        mvals.append(MASS[idx])
        
    #-- Printing values for corresponding index:       
        print("\nClosest Values to Selection:")
        print("\n\tRadius (km): {0}".format(R_val))    
        print("\tMass (M_sun): {0}".format(M_val))
        print("\tPressure (MeV/fm^3): {0}".format(P_val))
        print("\tCentral Energy Density (MeV/fm^3): {0}".format(E_val))
        print("\tk2: {0}".format(k2_val))
        print("\tCompactness: {0}".format(C_val))
        print("\tLambda: {0}".format(Lam_val))
        
        return(xx, yy)
        for sel in cursor.selections:
            cursor.remove_selection(sel)
        
    
    
    plt.tight_layout()
    plt.show()
    plt.close()

    #-- Outputting MR profile for selection (if a selection was made)
    if len(rvals) > 0:
        RADIUS = rvals[len(rvals)-1]
        MASS = mvals[len(mvals)-1]
    return(RADIUS, MASS)


#---Plotting MR Profile for Selection-------------------------------------------#   
def plotting2 (X, Y):    
    
    fig, ax = plt.subplots(1)
    graph = ax.plot(X, Y, 'c')
    ax.set_xlabel('Radius (km)', fontsize=15)
    ax.set_ylabel('Mass (M_sun)', fontsize=15)
    ax.set_title('MR Profile for Most Recent Selection\n(Click to Interact)', fontsize=15)
    ax.grid(True, which='major')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=15, length=10)
    ax.tick_params(axis='both', which='minor', length=5)

    cursor = mplcursors.cursor(graph)   # Allows user to interact w/ graph and get values
    @cursor.connect("add")
    def return_values(sel):
        tup = cursor.selections   # Gets x&y values for selection
        xx, yy = tup[0][1]
        return(xx, yy)
        for sel in cursor.selections:
            cursor.remove_selection(sel)

    plt.tight_layout()
    plt.show()
    plt.close()
    
#--Defining the first ODE (dP/dr) function--#
def f (x, y, z, s, args):
    k, c, e, dPdE = args
    A = z*c**2 + 4*pi*y*x**3
    num = -fact1 * (e+y) * A
    A = num / ( x**2 - 2*k*z*x )      #--A will have units J/m^4.
    return(A)

#--Defining the second ODE (dM/dr) function--#
def g (x, y, z, s, args):
    k, c, e, dPdE = args
    B = 4*pi*e*x**2/c**2              #--B will have units kg/m.
    return(B)  
    

#--Defining the third ODE (dyR/dr) function--#
def H (x, y, z, s, args):
    k, c, e, dPdE = args
    
    """#--conversions (refer the notes PDF, sec-1)--#
    x /= 1000                       #--in Km
    c /= 1000                       #--in Km/s
    k *= 1.476981739                #--in Km/M_sun
    z /= 1.98892*10**30             #--in M_sun 
    e *= 5.0278396266*10**(-28)
    y *= 5.0278396266*10**(-28)     #--in M_sun/km.s^2"""

    Fr = (x - 4*pi*(e-y)*(fact1)*x**3) / (x - 2*k*z)
    
    num1 = 4*pi*x * ( ( 5*e + 9*y + ((e+y)/dPdE) )*(fact1) - (6/(4*pi*x**2)) )
    den1 = x - 2*k*z
    part1 = num1/den1
    
    num2 = z*c**2 + 4*pi*y*x**3
    den2 = x**2 - 2*k*z*x
    part2 = 4*(fact1*num2/den2)**2
    Qr = part1 - part2
    
    H = (-1/x) * (s**2 + s*Fr + Qr*x**2) 
    
    return(H)

#--k2--#
def k2_value (yR, C):
    a = (8.0*C**5.0/5.0) * (1.0-2.0*C)**2.0 * (2.0+2.0*C*(yR-1.0)-yR)
    ba = 2.0*C * (6.0 - 3.0*yR + 3.0*C*(5.0*yR-8.0)) + 4.0*C**3.0 * (13.0 - 11.0*yR + C*(3.0*yR-2.0) + 2.0*C**2.0 * (1.0+yR))
    bb = 3.0*(1.0-2.0*C)**2.0 * (2.0 - yR + 2.0*C*(yR-1.0)) * log(1.0-2.0*C)
    return(a*(ba+bb)**(-1.0))

#---------------------------------------------------------------------------#    

#--Finding the matching value--#

def interpolateE (yi1, P, E): 

    for idx in range(1,len(P)-1):
        if yi1==P[idx]:
            e = E[idx]
            
        elif yi1>P[idx] and yi1<P[idx+1]:
            x0 = yi1
            x1, x2 = P[idx], P[idx+1]
            y1, y2 = E[idx], E[idx+1]
             
            e = ((y2-y1)*(x0-x1)/(x2-x1) ) + y1
            
    return(e)


def interpolateP (e0, P, E): 

    for idx in range(1,len(E)-1):
        if e0==E[idx]:
            p = P[idx]
            
        elif e0>E[idx] and e0<E[idx+1]:
            x0 = e0
            x1, x2 = E[idx], E[idx+1]
            y1, y2 = P[idx], P[idx+1]
             
            p = ((y2-y1)*(x0-x1)/(x2-x1) ) + y1
            
    return(p)
       
#------------------------------------------------------------------------#

#--RK4 Method--#
def rk4 (xi, yi, zi, si, h, args):

    k1 = f(xi, yi, zi, si, args)
    l1 = g(xi, yi, zi, si, args)
    j1 = H(xi, yi, zi, si, args)
    
    k2 = f(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    l2 = g(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    j2 = H(xi+h/2, yi+k1/2, zi+l1/2, si+j1/2, args)
    
    k3 = f(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    l3 = g(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    j3 = H(xi+h/2, yi+k2/2, zi+l2/2, si+j2/2, args)
    
    k4 = f(xi+h, yi+k3, zi+l3, si+j3, args)
    l4 = g(xi+h, yi+k3, zi+l3, si+j3, args)
    j4 = H(xi+h, yi+k3, zi+l3, si+j3, args)

    yi1 = yi + h/6*(k1 + 2*k2 + 2*k3 + k4)
    zi1 = zi + h/6*(l1 + 2*l2 + 2*l3 + l4)
    si1 = si + h/6*(j1 + 2*j2 + 2*j3 + j4)

    return(yi1, zi1, si1)
    

#--The TOV Solver Program--#
def eq_solve (L, e0, h, P0):
    
    P = L[:, 0]                          #--first column is pressure
    E = L[:, 1]                          #--second column is energy density
    
    fact = 1.60218*10**32                
    P, E = fact*P, fact*E                #--converted to SI units MeV/fm^3 to Pa=J/m^2. 

    dP = np.gradient(P)
    dE = np.gradient(E)
    dPdE = dP/dE                        #--this is the array dP/de.
    
    xMax = 30000  
    x0 = 0.1        #--in meter.

    
    m0 = (4/3)*pi*x0**3*e0/c**2           #--this is in Kg.
    
    p0 = interpolateP(e0, P, E)           #--these are also in desired units. 
    dPdE0 = interpolateP(e0, dPdE, E)     #--this gives the starting value of dPdE.
    s0 = 2

    P0.append(p0/fact)
    
    xi, yi, zi, si, ei, dPdEi = x0, p0, m0, s0, e0, dPdE0     #--initial conditions.

    mm = []
    rr = []
    
    for xi in frange(x0, xMax, h):
        args = [k, c, ei, dPdEi]
        yi, zi, si = rk4 (xi, yi, zi, si, h, args)
        rr.append(xi/1000)                  #-- units km
        mm.append(zi/(1.98847*10**30))      #--units M_sun
        
        if yi<0:
            break
        
        ei = interpolateE(yi, P, E)           #--thus e is already in the J/m^3.
        dPdEi = interpolateP(ei, dPdE, E)

    return(xi, yi, zi, si, mm, rr)
    return(P0)

#---------------------------------------------------------------------------#
    

def interpolate_vals(mass, radius, P, E, k2, C, Lam):
#-- Not yet compatible for mass and radius profiles (MASS & RADIUS) 
    mass_new = []
    radius_new = []
    P_new= []
    E_new = []
    k2_new = []
    C_new = []
    Lam_new = []
    
    for i in range(0,len(mass)-1):

        #-- Adding original val
        mass_new.append(mass[i])    
        radius_new.append(radius[i])
        P_new.append(P[i])
        E_new.append(E[i])
        k2_new.append(k2[i])
        C_new.append(C[i])
        Lam_new.append(Lam[i])

        #-- Adding interpolated val at midpoints (1/2^k)
        #--- k adjusts step-size
        k = 10;
        for j in range(1,2**k):
            mass_new.append(mass[i] + j*(mass[i+1]-mass[i])/2**k)
            radius_new.append(radius[i] + j*(radius[i+1]-radius[i])/2**k)
            P_new.append(P[i] + j*(P[i+1]-P[i])/2**k)
            E_new.append(E[i] + j*(E[i+1]-E[i])/2**k)
            k2_new.append(k2[i] + j*(k2[i+1]-k2[i])/2**k)
            C_new.append(C[i] + j*(C[i+1]-C[i])/2**k)
            Lam_new.append(Lam[i] + j*(Lam[i+1]-Lam[i])/2**k)


    #Reformatting
    mass = np.array(mass_new)
    radius = np.array(radius_new)
    P = np.array(P_new)
    E = np.array(E_new)
    k2 = np.array(k2_new)
    C = np.array(C_new)
    Lam = np.array(Lam_new)
    return(mass, radius, P, E, k2, C, Lam)



#---------------------------------------------------------------------------#        
#--Main Program--#
def main ():
    L = []
    while(len(L)==0 or len(L[0,:]) != 3):
        try:
            EoS_file = input('\nPlease give the name of the EoS file (eosaapr.dat):\t')
            while(EoS_file[len(EoS_file)-4:len(EoS_file)] != '.dat'):
                print('\n\tERROR -- Please enter .dat file\n')
                EoS_file = input('Please give the name of the EoS file (eosaapr.dat):\t')
            
            L = np.loadtxt(EoS_file)   #--Reading the EoS.
            if len(L[0,:]) !=3:
                print('\n\tERROR -- Incorrect dimensions\n')
 
        except Exception:
            print('\n\tERROR -- File name is not correct!\n')
            pass
        
        
    

    fact = 1.60218*10**32            #--to SI.
    e0 = float(input('Please give the value (eg: 153.27) of center energy density [MeV/fm^3]: '))    
    print("\n e0 in J/m^2 is: ", e0*fact, '\n')

    #h = float(input('Please give the radial step size (meters):\t'))
    h = 10                           #--radial step size in meter.

    
    N = int(input('Now give the number of stars (eg: 150): '))
    

    dE = 0
    radius, mass, yR, E0, P0 = [], [], [], [], []

    MASS = []           #-- Mass & Radius profiles (all values for each star)
    RADIUS = []

    #-- Solving ODE for each value of e0
    for i in range(0, N):
        try:
            e0 = e0 + dE
            E0.append(e0)
            X, Y, Z, S, mm, rr = eq_solve(L, e0*fact, h, P0)
            radius.append(X)
            mass.append(Z)
            yR.append(S)
            MASS.append(mm)
            RADIUS.append(rr)
            dE = 10
        
        except UnboundLocalError:
            print('\nThe max number of stars:', i, ' that can be plotted, is reached.')
            N = i
            break

   
    #--Defining compactness, k2, tidal polarizability (lambda)
    c = []
    k2 = []
    Lam = []
    for j in range(0, N):
        C = k*(mass[j]/radius[j])
        k2.append(k2_value(yR[j],C))
        Lam.append(2*(k2[j]*(radius[j]/(k*mass[j]))**5)/3)
        c.append(C)
    C = c;
        
    [radius, mass, yR, k2] = np.array(radius), np.array(mass), np.array(yR), np.array(k2)
    radius = radius/1000
    mass = mass/M_sun
    
   
   #--- Interpolation: (Not yet Compatible with MASS & RADIUS)
   # mass, radius, P0, E0, k2, C, Lam = interpolate_vals(mass, radius, P0, E0, k2, C, Lam)

 

    print('\n For the given EoS, we have:\n  Maximum Mass -> ',
          max(mass), ' M_sun', '\n  Corresponding Radius  -> ',
          float(radius[np.where(mass==max(mass))]), ' Km', 
          '\n  Corresponding k2 -> ', float(k2[np.where(mass==max(mass))]), '\n')   

#-------------------------Saving Output--------------------------#
    ans1 = input('Do you want to store the output? (Y/N):\t')
    while (ans1 !='Y' and ans1 !='N'):
        print('\n\tERROR -- Please choose Y or N\n')
        ans1 = input('Do you want to store the output? (Y/N):\t')
        
    if ans1=='Y':
        fname1 = input('Give the table name in which to store the output(.txt):\t')
        while(fname1[len(fname1)-4:len(fname1)] != '.txt'):
            print('\n\tERROR -- Please save as .txt file\n')
            fname1 = input('Give the table name in which to store the output(.txt):\t')

        Matrix = np.zeros((len(mass),7))
        for i in range(0,len(mass)):
            Matrix[i][0] = radius[i]
            Matrix[i][1] = mass[i]
            Matrix[i][2] = k2[i]
            Matrix[i][3] = C[i]
            Matrix[i][4] = P0[i]
            Matrix[i][5] = E0[i]
            Matrix[i][6] = Lam[i]

        np.savetxt(fname1, Matrix)
        print('\n\tSaved as matrix {0}x7 with columns: R,M,k2,C,P,E,Lam \n'.format(len(mass)))                      
 #---------------------------------------------------------------------##              
        
    ans2 = input('\nDo you want a plot? (Y/N):\t')
    while (ans2 !='Y' and ans2 !='N'):
        print('\n\tERROR -- Please choose Y or N\n')
        ans2 = input('Do you want a plot? (Y/N):\t')
    
    while ans2=='Y':
        a1 = ['radius', 'mass', 'k2', 'c', 'P', 'E', 'lambda']
        a2 = ['line', 'dot', 'dash']

        query1 = input('On x-axis (radius/mass/k2/c/P/E/lambda):\t')
        while query1 not in a1:
            print('\n\tERROR -- Please choose from the options listed\n')
            query1 = input('On x-axis (radius/mass/k2/c/P/E/lambda):\t')
            
        query2 = input('On y-axis (radius/mass/k2/c/P/E/lambda):\t')
        while query2 not in a1:
            print('\n\tERROR -- Please choose from the options listed\n')
            query2 = input('On y-axis (radius/mass/k2/c/P/E/lambda):\t')
            
        query3 = input('Which marker (line/dot/dash):\t')
        while query3 not in a2:
            print('\n\tERROR -- Please choose from the options listed\n')
            query3 = input('Which marker (line/dot/dash):\t')
        
        title = query2 + '-' + query1 + ' Plot' + ' (Click to Interact)'
        
        if query1=='radius':
            x = radius
            xname = 'Radius (Km)'
        if query2=='radius':
            y = radius
            yname = 'Radius (Km)'
        if query1=='mass':
            x = mass
            xname = 'Mass (M_sun)'
        if query1 == 'c':
            x = C
            xname = "Compactness"
        if query1 == 'P':
            x = P0
            xname = 'Pressure (MeV/fm^3)'
        if query1 == 'E':
            x = E0
            xname = 'Energy Density (MeV/fm^3)'
        if query1 == 'lambda':
            x = Lam
            xname = 'Lambda'
    
        if query2=='mass':
            y = mass
            yname = 'Mass (M_sun)'
        if query1=='k2':
            x = k2
            xname = 'k2'
        if query2=='k2':
            y = k2
            yname = 'k2'
        if query2=='c':
            y = C
            yname = 'Compactness'
        if query2 == 'P':
            y = P0
            yname = 'Pressure (MeV/fm^3)'
        if query2 == 'E':
            y = E0
            yname = 'Energy Density (MeV/fm^3)'
        if query2 == 'lambda':
            y = Lam
            yname = 'Lambda'
    
        
        if query3=='line':
            marker = 'c'
        if query3=='dot':
            marker = '.'
        if query3=='dash':
            marker = '--'    

        RProfile, MProfile = [], []
        RProfile, MProfile = plotting(x, y, mass, radius, P0, E0, k2, C, Lam, RADIUS, MASS, marker, xname, yname, title)

        if len(RProfile) > N:       #-- Only excecutes if a selection was made from previous plot
            plotting2 (RProfile, MProfile)
       
        ans2 = input('\nDo you want another plot? (Y/N):\t')
        while (ans2 !='Y' and ans2 !='N'):
            print('\n\tERROR -- Please choose Y or N\n')
            ans2 = input('Do you want another plot? (Y/N):\t')
    
    answer = str(input('Do you want to use another EoS? (Y/N):\t'))
    while (answer !='Y' and answer !='N'):
            print('\n\tERROR -- Please choose Y or N\n')
            answer = input('Do you want to use another EoS? (Y/N):\t')
            
    while answer=='Y':
        answer='N'
        main()


main()
        

   
#####################--End of Program--##################################
#########################################################################
