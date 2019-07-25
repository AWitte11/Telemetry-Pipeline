"""NollEstimate.py - Estimation of Fried's Parameter and seeing angle given DM commands

Language: Python 3

Written by: Alex Witte for the 2019 University of Hawaii REU Program

Mentor: Dr. Charlotte Bond

Last Updated:
    07/24/2019
"""


import math
import numpy as np
from scipy import io
from matplotlib import pyplot as plt
from astropy.io import fits
import statistics as stat
import Decompose as dc
from zernike import zernikeArray
import sys
 


def fried(data,jmin,jmax,D,nz,rad,save_dir,save_plot,sim,plot = True):
    """Estimates the Fried Parameter, using Noll's approximation, from the variances of the Zernike coefficients
    
    Arguments:
        data(numpy array): Array of raw data
        jmin(int): minimum Noll index used in fitting the theoretical curve
        jmax(int): maximum Noll index used in fitting the theoretical curve
        D(int): diameter of the wavefront being decomposed
        nz(int): number of Zernike Modes to use
        rad(int): Diameter of inner circle to decompose in actuators
        save_dir(str): String containing the directory to save the plot and data
        save_plot(Bool): Whether or not to save plot
        plot(Bool): Whether or not to plot data
        
    
    Returns:
        see_arc(float): An estimate of the seeing angle in arcseconds derived from the Zernike decomposition
        r0(float):  An estimate of the Fried parameter in arcseconds derived from the Zernike decomposition
    """
    if sim == True:
        outer = convert_sim(data)
        print("Conversion Successful")
        inner = inner_circle(outer,rad)
        print("Inner Circle Success")
    
    elif sim == False:
    #        Applies influence funtion
        outer = influence(data)
        rad = 5 * rad
        print('Influence Success')
        
        
        ###cuts out the desired inner circle for Zernike decomposition
        inner = inner_circle(outer,rad)    
        print("Inner Circle Success")
            
    
    
#   Decompose all the data        
    master,difs = dc.decomp_all(inner,nz,rad)
    
    #Set up to find variance
    var_list = []
    var_last = 0
    
#   Find variances
    for row in master:
        mean = stat.mean(row)
        s = 0
        for element in row:
            s += (element - mean)**2
        var_last = (s/len(row))
        var_list.append(var_last)        
        
    print("Done with variances") 
#    Calculation of error in the variances
    var_err = stat.variance(difs)
    
#   Coefficients of Noll's approximation for J <= 10
    an = [1.030,0.582,0.134,0.111,0.088,0.065,0.059,0.053,0.046,0.04]


    lastdif = 10**20
#    Scan through r0's untill the correct fit is found
    for r in np.linspace(.001,0.2,20000):
        noll_fit = []
        
#        Create the theoretical curves according to Noll's approximations 
        for j in range(0,10):
            noll_fit.append(an[j]*((D/r)**(5/3)))               
        for J in range(11,nz+1):
            noll_fit.append((0.2944*(J**(math.sqrt(3)/-2))*((D/r)**(5/3))))
            
            
#       Take the differential of the curves  
        fit_diff = np.diff(noll_fit)
        
        
#       Calculate the difference between the data and the fit line
        seperations = []
        for j in range(6,75):
            seperations.append(abs(abs(fit_diff[j]) - var_list[j]))

        if np.mean(seperations) > lastdif:
            r0 = r
            break
        
        lastdif = np.mean(seperations)
    
    
    aprox = []
    
    for j in range(0,10):
            aprox.append(an[j]*((D/r0)**(5/3)))               
    for J in range(11,nz+1):
            aprox.append((0.2944*(J**(math.sqrt(3)/-2))*((D/r0)**(5/3))))
            
            
    aprox_difs = np.diff(aprox)
    aprox_plot = []
    for element in aprox_difs:
        aprox_plot.append(abs(element))
    
#   Calculate error on r0
    r_err = (0.2944*(J**(math.sqrt(3)/-2))*((D/var_err)**(5/3))) 
   
    
#    Convert r0(m) to seeing angle(arcsecs) @ 500nm
    see_radn = (.98)*500e-9*(1/r0)
    see_arcn = see_radn * (3600*180)/(math.pi)
    
#    Print results
    print('r n: ',format(r0,'.7g'))
    print('r err n: ',format(r_err,'.7g'))
    print('see n: ',format(see_arcn,'.7g'))
    
#   Plot Results
    if plot == True:           
        plt.yscale('log')        
        plt.plot(var_list,'b')
        plt.plot(aprox_plot,'k',label = 'Fit r0 = '+format(r0,'.5g')+' m')
        plt.title('Variance vs. Zernike Mode On-Sky') 
        plt.legend()
        plt.xlabel('Zernike Mode')
        plt.ylabel('Variance (rad^2)')
        
#        Save Results
        if save_plot == True:
            plt.savefig(save_dir.strip('"')+'//variance_vs_Mode.jpeg')
            
        plt.legend()
        plt.show()

#   Function returns seeing estimmate as well as r0 estimate
    return see_arcn,r0



def inner_circle(data,diam):
    """Pulls inner circle of actuators for Zernike decomp
    
    Arguments:
        data(numpy array): raw phase map
        diam(int): Diameter of the inner circle, in actuators, to extract
    
    """
    
    fin = []
    
#    Generates the control array with the desired shape
    control = np.ndarray.flatten(zernikeArray(1,diam))
    
    
#    Set edges of the inner circle
    minimum = ((len(data[0])-diam)//2)
    maximum = (((len(data[0])-diam)//2)+diam)
    
    
    for frame in data:
        hold = np.zeros((diam,diam))
        
#        Pulls the inner radius from the raw data
        inner = np.ndarray.flatten(frame[minimum:maximum,minimum:maximum])
        
        
#        Conforms the inner circle to the control shape
        n = 0
        for j in range(0,diam):
            for i in range(0,diam):                
                if control[n] != 0:
                    hold[j][i] = inner[n]                    
                n+=1    
        fin.append(hold)
          
    return fin


def convert_sim(filename):
    """Converts simulated data into the same fromat as the on sky data
    
    Arguments:
        filename(string): file containing the simulated data
        
    Returns:
        final(numpy array): Array containing the simulated data in the on sky format
    """
    
    with fits.open(filename) as hdul:
            data = hdul[0].data
            
    final = []
    

#    for each frame
    for frame in data:
        i = 0
       
        hold = np.zeros((21,21))

#       Converts to the same shape as the on-sky as well as converting from volts to phase
        for row in range(0,21):
            for colm in range(0,21):
                hold[row][colm] = 2*math.pi*(1/500)*(10e8*frame[i])
                i += 1

        final.append(hold)
        
    
    return final



def influence(data):
    """Applies the influence function to the DM commands
    
    Arguments: 
        data(numpy Array): raw DM commands
        
    Returns:
        final(numpy array): True shape of the wavefront
    """
    
#    Load in pre-generated Influence function
    influence = io.loadmat("Influence 101.mat")['ans']
    
#    Set indices of DM commands in larger shape array
    indexes = np.linspace(0,100,21)
    z = [ (a,b) for a in indexes for b in indexes ]



    final = []
    framenum = 1

#   Takes each frame of 21x21 DM commands and scales it by a factor of 5 to allow for the shape of the influence function       
    for frame in data:  
        data_big = np.zeros_like(influence[:,:,1])
        for row in range(0,len(frame)):
            for column in range(0,len(frame)):
                data_big[int(5 * row)][int(5 * column)] = frame[row][column] 
                
        
        hold = np.zeros_like(data_big)
        
#        For each of the acutators, apply the influence function
#        Sum these influences together to get the true shape of the DM
        for mode in range(0,441):
            column,row = z[mode] 
            weight = data_big[int(row),int(column)]
            weighted = np.multiply(weight,influence[:,:,mode])
            hold = np.add(hold,weighted)
        
#        Convert from volts to nm
        final.append(np.multiply(2*math.pi*(480/500),hold))
        
        if framenum == (len(data)//2):
            print("Halfway of Influence")
        framenum+=1
        
        
    return final
