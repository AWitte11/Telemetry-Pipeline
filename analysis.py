"""Analysis.py - Several analysis methods

Language: Python 3

Written by: Alex Witte for the 2019 University of Hawaii REU Program

Mentor: Dr. Charlotte Bond

Last Updated:
    07/24/2019
"""

import math
import numpy as np
from matplotlib import pyplot as plt
from scipy import io
import Decompose as dc
import NollEstimate as ne
import statistics as stat



def rms(raw,save_dir,title,Plot=True):
    """Plots the Root Mean Square of the DM/residual WF for each frame
    
    Arguments:
        raw(numpy array): raw data
        save_dir(str): String containing the directory to save the plot and data
        title(string): name of type of file being analized (Residual WF vs DMC)
        Plot(bool): Plot or not
        
    Returns:
        final(list): list of RMS values for each frame
    """
    
#    Convert to nm
    data = np.multiply(raw,600)
 
    final = []
    print('Generating RMS')
#    Load in the WFS pupil to act as a guide. If the pixle isn't illuminated on the pupil, the point is ignored
    control_array = np.load('WFSPupil.npy')


#    For each frame in the data: 
#        1)Check if the point is on the pupil
#        3)Use the valid points to create the rms for the frame
    for point in data:
    
        j = 0 #row
        i = 0 #column
        val = 0
        numb = 0
        for j in range(0,len(point[0])):
            for i in range(0,len(point[0])):           
                if control_array[j,i] != 0:
                    val += (point[j,i])**2
                    numb += 1
        final.append(math.sqrt(((val))/numb))

    print("Begin Plot")

#    Create and save the plots for the RMS
    if 'dmc' or 'DMC' in title:
        plt.plot(final,zorder = 0)
        plt.title('DMC RMS') 
        plt.xlabel('Frame')
        plt.ylabel('RMS (nm)')
        plt.hlines(np.average(final),0,len(final),'r',label = 'Avg: '+format(np.average(final),'.5g')+ ' nm',zorder = 10)
        plt.legend()
        plt.savefig(save_dir+'\\DMC_RMS.pdf')
        
        
    elif 'residual' in title:
        plt.plot(final,zorder = 0)
        plt.title('Residual RMS') 
        plt.xlabel('Frame')
        plt.ylabel('RMS (nm)')
        plt.hlines(np.average(final),0,len(final),'r',label = 'Avg: '+format(np.average(final),'.5g')+ ' nm',zorder = 10)
        plt.legend()
        plt.savefig(save_dir+'\\residual_RMS.pdf')
        
    if Plot == True:
        plt.show()
        
    
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



def low_order_analysis(data,rad,nz):
    """Collects data on the low order Zernike coefficients
    
    Arguments:
        data(numpy array): raw data
        
    """
    
    #        Applies influence funtion
    outer = influence(data)
    rad = 5 * rad
    print('Influence Success')
    
    
    ###cuts out the desired inner circle for Zernike decomposition
    inner = ne.inner_circle(outer,rad)    
    print("Inner Circle Success")
            
    
    
#   Decompose all the data        
    master,difs = dc.decomp_all(inner,nz,rad)
    
    #    Print stats on coefficients
    for j in range(0,11):
        print('J= ',j)
        print('Mean coefficient: ',stat.mean(master[j]))
        print('Std. Deviation: ',stat.stdev(master[j]))        
