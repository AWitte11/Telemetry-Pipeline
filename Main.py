"""Main.py - Main run file for telemetry analysis

Language: Python 3

Written by: Alex Witte for the 2019 University of Hawaii REU Program

Mentor: Dr. Charlotte Bond

Last Updated:
    07/24/2019
"""


import NollEstimate as NE
import analysis
import numpy as np
from astropy.io import fits


###Parameters###         
###Min and max modes to approximate over
jmin = 6
jmax = 75

###physical diameter of pupil (if using the inner 17x17 actuators D = 9m)
D = 9

###Number of Zernike Modes to produce 
nz = 100

###Diameter of inner circle to decompose in actuators 
diam = 17



if __name__ == '__main__':
    
    issim = input("Is this simulated data? Y/N ")
              

#Select file for analysis and open it
    filename = input("File path: ")     
    if filename[0] == '"':
        filename = filename.strip('"') 
    
    if issim == 'Y':
        with fits.open(filename) as hdul:
            data = hdul[0].data
    elif issim == 'N':
        data = np.load(filename)

    print('Load Successful')
    
    
    
#Set up save directory
    save_dir = input("Save Directory (if none = False): ")    
    if save_dir[0] == '"':
        save_dir = save_dir.strip('"')
        
    save = True    
    if save_dir == "False":
        save = False
     
      
#        Decide the analysis to do
    title = input("DMC or residual?: ")        
    rms = input("Calculate RMS? (Y/N): ")    
    fri = input("Estimate seeing/r0? (Y/N): ")
    low = input("Analize lower order modes? (Y/N): ")
    
    
   #Preform analysis 
    if rms == 'Y':
        analysis.rms(data,save_dir,title,Plot=True)
        
    if fri == 'Y':
        
       
        if issim == 'Y':
            NE.fried(data,jmin,jmax,D,nz,diam,save_dir,save,True,plot = True)
        elif issim == 'N':
            NE.fried(data,jmin,jmax,D,nz,diam,save_dir,save,False,plot = True)

    if low =='Y':
        analysis.low_order_analysis(data,diam,nz)