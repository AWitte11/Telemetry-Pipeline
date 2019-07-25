# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:57:09 2019

@author: alexw
"""

import numpy as np
from zernike import zernikeArray
from matplotlib import pyplot as plt


def zern_array(nz,rad):    
    """Creates the inverse matrix used in zern decomposition
    
    Returns:
        final(numpy array): Inverse Matrix of Zernike basis
    """         
    
    #Create Zernike polynomials
    master_zern = zernikeArray(nz,rad)
    zern_matrix = np.zeros(((rad**2),len(master_zern)))
    
    
    j = 0 #row
    i = 0 #column
    
#    Assembles a matrix in which each column corresponds to a Zernike polynomial
#    and each row corresponds to a point in the DM command/residual WF
    for mode in master_zern:
        flat = mode.flatten()  
        for j in range(0,(rad**2)):
            zern_matrix[j][i] = flat[j]
        i+=1
     
    #create inverse matrix
    final = np.linalg.pinv(zern_matrix,rcond = 1e-10)
    
    return final



def zern_decomp(data,conversion, Plot = False):
    """Decomposes the frame into Zernike basis
    
    Arguments:
        data(numpy array): frame of data
        conversion(numpy array): Inverse Matrix computed by zern_array()
        
    Returns:
        coe(numpy vector): Vector containing the coefficients for the Zernike Ploynomials
    """
   
    #Flattens data matrix to a vctor
    flat = np.matrix.flatten(data)
    
    #Compute the coefficients by matrix multiplication
    coe = np.matmul(conversion,flat) 
    
#    Plot the decomposition of the frame
    if Plot == True:        
        plt.bar(range(0,len(coe)),coe)
        plt.title('Zernike Decomposition') 
        plt.xlabel('Zernike Mode')
        plt.ylabel('Zernike Coefficient')
        plt.show()
    
    return coe



def decomp_all(data,nz,rad):
    """ Decomposes all frames into Zernike polynomials
    
    Arguments:
        data(numpy array): raw array containing each frame of the sensor reading
        
    Returns:
        master_array(numpy array): Each column is a different frame with each row corresponding to a
        a different polynomial coefficient
    """
    
    master_array = np.zeros((nz,len(data)))
    conversion = zern_array(nz,rad)
    
    point_dif = []
#    find the coefficients for each frame in the raw data
    for index in range(0,len(data)):
        point = data[index]
        coef = zern_decomp(point,conversion)
        
        reconstructed = np.dot(coef,conversion) 
            
        dif = []
        x = 0   
        
    #    Finds the difference between the reconstructed and the actual wavefront
        for row in range(0,len(point)):
            for col in range(0,len(point[row])):
                dif.append((reconstructed[x]-point[row][col]))
                x+=1
        
        point_dif.append(np.mean(dif))        
        
#        arrange the coeficients in columns of a master array
        for x in range(0,len(coef)):
            master_array[x][index] = coef[x]
            
    return master_array,point_dif

