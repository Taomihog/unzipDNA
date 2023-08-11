# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:15:18 2020

@author: xgao
"""
import numpy as np
#Units are pN.nm.s if not specified.

#constants
Kelvin = 298.15
Joul = 4184
Avogadro = 6.022E+23
Boltzmann=0.0138065
#kT=Kelvin*Boltzmann

#DNA parameters
#LPDS = 51.97
#KDS = 1318
#L0DS = 0.338
#LPSS = 0.765
#KSS = 470
#L0SS = 0.554
    
def getEDNA(seq,salt_mM,kT):
#    print (len(seq))
    LUTseq = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    LUTdH = [-7.28, -5.8, -5.21, -4.63, -8.96, -8.57, -9.66, -5.21, -8.16, -10.1, -8.57, -5.8, -8.31, -8.16, -8.96, -7.28]
    LUTdS = [-20.28, -14.46, -12.89, -11.62, -24.48, -22.3, -24.43, -12.89, -22.46, -25.96, -22.3, -14.46, -25.06, -22.46, -24.48, -20.28]
    LUTm = [0.145, 0.099, 0.07, 0.117, 0.091, 0.063, 0.132, 0.07, 0.155, 0.079, 0.063, 0.099, 0.091, 0.155, 0.091, 0.145]
    index=np.zeros(len(seq)-1,dtype = np.int8)
    energy=np.zeros(len(seq),dtype = np.float)
    eff_salt = np.log(salt_mM*0.001)/298 #it is not equal to the room temperature
    for i in range(len(index)):
        for j,s in enumerate(LUTseq):
            if LUTseq[j]==seq[i:i+2]:
                break
#            print ('unknown nucleotide!!')
        index[i]=j
        energy[i+1]=energy[i]+LUTdH[j]-(LUTm[j]*eff_salt+LUTdS[j]*0.001)*Kelvin
    energy=-kT*energy*1e21*Joul/Avogadro/Boltzmann/Kelvin
    for i in range (10) :
        print (energy[i])
    print ("...")
    print (energy[-1])
#    np.savetxt('energy.txt',energy.T,delimiter = '\t',header ='energy')
    return energy

def Z_lev (kpillar,F):
    return F/kpillar
def Z_lev_df (kpillar,F):
    return 1/kpillar

def phi_ss (F,LPSS,KSS,kT):
    return (1/np.tanh(2*F*LPSS/kT)-kT/(2*F*LPSS))*(1+F/KSS)
def phi_ss_df (F,LPSS,KSS,kT):
    a = np.tanh(2*F*LPSS/kT)
    return (1/a-kT/(2*F*LPSS))/KSS+((1-1/(a*a))*2*LPSS/kT+kT/(2*F*F*LPSS))*(1+F/KSS)

#its better to use WLC here? KDS may not be necessary but I am not sure
def phi_ds (F,LPDS,KDS,kT):
    return 1-0.5*np.sqrt(kT/(F*LPDS))+F/KDS
def phi_ds_df (F,LPDS,KDS,kT):
    return 1+0.25*np.sqrt(kT/(F*F*F*LPDS))+1/KDS


#def Zlinker(F,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
#    return Z_lev(kpillar,F)+2*jss*L0SS*phi_ss(F,LPSS,KSS,kT)+jds*L0DS*phi_ds(F,LPDS,KDS,kT)#There is a factor of 2 in ssDNA length!!!
#def Jac_delta_Zlinker(F,z0,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
#    return -(Z_lev_df(kpillar,F)+2*jss*phi_ss_df(F,LPSS,KSS,kT)+jds*L0DS*phi_ds_df (F,LPDS,KDS,kT)) 

#def Zlinker(F,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
#    return F/kpillar+2*jss*L0SS*((1/np.tanh(2*F*LPSS/kT)-kT/(2*F*LPSS))*(1+F/KSS))+\
#            jds*L0DS*(1-0.5*np.sqrt(kT/(F*LPDS))+F/KDS)
#            #There is a factor of 2 in ssDNA length!!!
def Zlinker(F,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
    return F/DNA_parameters[0]+\
                         2*DNA_parameters[4]*DNA_parameters[3]*((1/np.tanh(2*F*DNA_parameters[1]/DNA_parameters[9])-DNA_parameters[9]/(2*F*DNA_parameters[1]))*(1+F/DNA_parameters[2]))+\
                         DNA_parameters[8]*DNA_parameters[7]*(1-0.5*np.sqrt(DNA_parameters[9]/(F*DNA_parameters[5]))+F/DNA_parameters[6])
            #There is a factor of 2 in ssDNA length!!!

def Jac_delta_Zlinker(F,z0,DNA_parameters):
    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
    a = np.tanh(2*F*LPSS/kT)
    return -(1/kpillar+\
            2*jss*((1/a-kT/(2*F*LPSS))/KSS+((1-1/(a*a))*2*LPSS/kT+kT/(2*F*F*LPSS))*(1+F/KSS))+\
            jds*L0DS*(1+0.25*np.sqrt(kT/(F*F*F*LPDS))+1/KDS)) 

#For test, the result is correct
def plotForce(DNA_parameters):
    N = 700
    force = np.empty(N)
    extension = np.empty(N)
    for F in range(N):
        force[F]=F*0.1+0.1
        extension[F] = Zlinker(F*0.1+0.1,DNA_parameters)
    np.savetxt('test.txt',np.vstack((extension, force)).T,delimiter = '\t',header ='ext\tforce')
    import matplotlib.pyplot as plt
    plt.plot(extension,force,'o-')
    plt.show()


#def delta_Zlinker(F,z0,DNA_parameters):
#    return z0-Zlinker(F,DNA_parameters)
def delta_Zlinker(F,z0,DNA_parameters):
    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
    return z0-(F/kpillar+2*jss*L0SS*((1/np.tanh(2*F*LPSS/kT)-kT/(2*F*LPSS))*(1+F/KSS))+\
            jds*L0DS*(1-0.5*np.sqrt(kT/(F*LPDS))+F/KDS))

#def delta_Zlinker2(F,z0,DNA_parameters):
#    return np.square(z0-Zlinker(F,DNA_parameters))
#def delta_Zlinker2(F,z0,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
#    #1      2    3   4    5   6    7   8    9   10
#    return np.square(z0-(F/kpillar+2*jss*L0SS*((1/np.tanh(2*F*LPSS/kT)-kT/(2*F*LPSS))*(1+F/KSS))+\
#            jds*L0DS*(1-0.5*np.sqrt(kT/(F*LPDS))+F/KDS)))
def delta_Zlinker2(F,z0,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
    #0      1    2   3    4   5    6   7    8   9    DNA_parameters[9]
    return np.square(z0-(F/DNA_parameters[0]+\
                         2*DNA_parameters[4]*DNA_parameters[3]*((1/np.tanh(2*F*DNA_parameters[1]/DNA_parameters[9])-DNA_parameters[9]/(2*F*DNA_parameters[1]))*(1+F/DNA_parameters[2]))+\
                         DNA_parameters[8]*DNA_parameters[7]*(1-0.5*np.sqrt(DNA_parameters[9]/(F*DNA_parameters[5]))+F/DNA_parameters[6])))

from scipy.optimize import minimize
#from scipy.optimize import root_scalar #Do not know why there is no such module

def FindForce(z0,DNA_parameters):
#    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
    fguess=10
    #This minimization can be further optimized because the slope can be exactily calculated, Newton method should be better
    
#    res = optimize.root_scalar(delta_Zlinker, args = (z0,DNA_parameters), fprime=Jac_delta_Zlinker, method='newton',x0 = fguess)
##    print (res)
#    return res.root
    
    res = minimize(delta_Zlinker2, fguess, args = (z0,DNA_parameters), method='Nelder-Mead', tol=1e-3)#Nelder-Mead
#    print(res)
    return res.x[0]

from scipy.integrate import quad
#def Zlinker_integ(F,kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT):
#    DNA_parameters = kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT
#    return Zlinker(F,DNA_parameters)
def Zlinker_integ(F,kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT):
#    DNA_parameters = kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT
    return Zlinker(F,(kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT))

'''has been integrated into getEtot
def getElinker (z0,DNA_parameters):
    f = FindForce(z0,DNA_parameters)
#    print (f)
#    return f*z0-quad(Zlinker_integ, 0.1, f, args=(DNA_parameters))[0]
    return f,f*z0-quad(Zlinker_integ, 0.1, f, args=(DNA_parameters))[0]
'''
def getEtot (z0,DNA_parameters,EDNA):
    #DNA_parameters = kpillar,LPSS,KSS,L0SS,j,LPDS,KDS,L0DS,jds,kT
    if DNA_parameters[4] <0: 
        DNA_parameters= DNA_parameters[:4] + (0,) + DNA_parameters[5:]
        f = FindForce(z0,DNA_parameters)
        Elinker = f*z0-quad(Zlinker_integ, 0.1, f, args=(DNA_parameters))[0]
        return f, Elinker, Elinker, 0.0
    #return force, Etot, Elinker, EDNA
    else:
        f = FindForce(z0,DNA_parameters)
        Elinker = f*z0-quad(Zlinker_integ, 0.1, f, args=(DNA_parameters))[0]
        return f, EDNA[DNA_parameters[4]]+Elinker, Elinker, EDNA[DNA_parameters[4]]

def E_ss (F,lp,K,l0):
    #Smith 1996 science (?need ref.)
    return 0

def E_ds (F,lp,K,l0):
    #Odijk 1995 macromolecules
    return 0


def calculate (EDNA, DNA_parameters):
    j_max =len(EDNA)-1
    print ('j_max = ', j_max)
    kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT = DNA_parameters
#    plotForce(DNA_parameters)

    #start calculation
    z0_init = jds*L0DS+10/kpillar #10 pN as the start point
    z_res = 5.0 #nm
    maxstep = 1000#1711
    
    searchLength = 500
    minEnergyDiff = 13*kT #if larger than this, stop searching
    jmin = 0
    
    print ('Starting at extension (nm) and force (pN):',z0_init, FindForce(z0_init,DNA_parameters))
    zzjjff = np.zeros((maxstep,7))
    
    j = 0
    for index_z in range (maxstep):

        TotSteps = 0
        P = 0
        FP = 0
        FFP = 0
        FFP = 0
        jP = 0
        jjP = 0
        jjP = 0
        
        
        #find j of the minimal energy
        z0 = z0_init+index_z*z_res
        Etot_prev=np.inf
        dEtot = -1
        j-=1
        while (dEtot<=0):
            j+=1
            if j >= j_max: break
            DNA_parameters = kpillar,LPSS,KSS,L0SS,j,LPDS,KDS,L0DS,jds,kT
            f, Etot, Elinker, EDNAj = getEtot (z0,DNA_parameters,EDNA)
            TotSteps+=1
            dEtot = Etot-Etot_prev#negative if Etot decrease
            Etot_prev = Etot
        dEtot = -1
        j+=1
        if j >= j_max: break
        while (dEtot<=0):
            j-=1
            if j == 0: break
            DNA_parameters = kpillar,LPSS,KSS,L0SS,j,LPDS,KDS,L0DS,jds,kT
            f, Etot, Elinker, EDNAj = getEtot (z0,DNA_parameters,EDNA)
            TotSteps+=1
            dEtot = Etot-Etot_prev#negative if Etot decrease
            Etot_prev = Etot
        
        j+=1
        if j >= j_max: break
        DNA_parameters = kpillar,LPSS,KSS,L0SS,j,LPDS,KDS,L0DS,jds,kT
        f, Etot, Elinker, EDNAj = getEtot (z0,DNA_parameters,EDNA)
        TotSteps+=1
        
        Etot_min = Etot
        P = 1
        FP = f
        FFP = f*f
        jP = j
        jjP = j*j
        
        Etot_new_min = Etot_min
        for i in range(searchLength):#calculate a bunch of positive j
            i+=1
            if j+i >j_max: break
            DNA_parameters = kpillar,LPSS,KSS,L0SS,j+i,LPDS,KDS,L0DS,jds,kT
            f, Etot, Elinker, EDNAj = getEtot (z0,DNA_parameters,EDNA)
            TotSteps+=1
            Ptemp = np.exp(-(Etot-Etot_min)/kT)
            P += Ptemp
            FP += f*Ptemp
            FFP += f*f*Ptemp
            jP += (j+i)*Ptemp
            jjP += (j+i)*(j+i)*Ptemp
            if Etot-Etot_min > minEnergyDiff: break
            if Etot<Etot_new_min: 
                Etot_new_min = Etot
                jmin = j+i
            
        for i in range(searchLength):#then calculate a bunch of negative j
            i= 0-i-1
            if (j+i) <0: break
            DNA_parameters = kpillar,LPSS,KSS,L0SS,j+i,LPDS,KDS,L0DS,jds,kT
            f, Etot, Elinker, EDNAj = getEtot (z0,DNA_parameters,EDNA)
            TotSteps+=1
            Ptemp = np.exp(-(Etot-Etot_min)/kT)
            P += Ptemp
            FP += f*Ptemp
            FFP += f*f*Ptemp
            jP += (j+i)*Ptemp
            jjP += (j+i)*(j+i)*Ptemp
            if Etot-Etot_min > minEnergyDiff: break
            if Etot<Etot_new_min: 
                Etot_new_min = Etot
                jmin = j+i
        j = jmin
        print ('\rindex_z = ',index_z,'j_min = ',j,'Total steps',TotSteps)
        
        FP /= P
        FFP /= P
        FFP -=FP*FP
        jP /= P
        jjP /= P
        jjP -= jP*jP

        zzjjff[index_z,:]=[z0-FP/kpillar,FFP/kpillar,jP,jjP,FP,FFP,TotSteps]
    zzjjff = zzjjff[:index_z,:]#if break, the last element will not be inserted
    return zzjjff

if __name__ == "__main__":   
    with open('exampleData.txt', 'r') as file:
        sequence = file.read().replace('\n', '')
    
    kT = Kelvin*Boltzmann
    salt_mM =100
    energy = getEDNA(sequence,salt_mM,kT)
    
    #DNA_parameters and kT
    LPDS = 51.97 #dsDNA persistence length
    KDS = 1318 #dsDNA stretch modulus
    L0DS = 0.338 #dsDNA contour length per bp
    LPSS = 0.765 #ssDNA persistence length
    KSS = 470 #ssDNA stretch modulus
    L0SS = 0.554 #ssDNA contour length per nt
    jss = 0 #ssDNA #nt
    jds = 2200 #dsDNA #bp
    kpillar = 0.07406 #z stiffness of AOT (Here, power = 200 mW)
    parameters = kpillar,LPSS,KSS,L0SS,jss,LPDS,KDS,L0DS,jds,kT
    
    File_Name_Result = 'result.txt'
    np.savetxt(File_Name_Result,calculate(energy, parameters),delimiter = '\t',\
               header ='zMean\tzSD\tjMean\tjSD\tFMean\tFSD\tTotSteps')