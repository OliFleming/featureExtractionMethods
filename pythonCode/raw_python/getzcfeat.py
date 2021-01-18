import math
import numpy as np
from scipy.signal import lfilter

# GETZCFEAT - Gets the ZC feature (zero crossing).
#
# Syntax:  [t,s] = WaveformGen(time,interval,amp,frequency,phase,dc,noise,type)
#
# feat = getsscfeat(x,deadzone,winsize,wininc,datawin,dispstatus)
#
# Matlab Original Author Adrian Chan
#
# This function computes the ZC feature of the signals in x,
# which are stored in columns.
#
# The signals in x are divided into multiple windows of size
# winsize and the windows are space wininc apart.
#
# Inputs
#    x: 		columns of signals
#    deadzone:  +/- zone signal must cross to be considered a deadzone
#    winsize:	window size (length of x)
#    wininc:	spacing of the windows (winsize)
#    datawin:   window for data (e.g. Hamming, default rectangular)
#               must have dimensions of (winsize,1)
#    dispstatus:zero for no waitbar (default)
#
# Outputs
#    feat:     RMS value in a 2 dimensional matrix
#              dim1 window
#              dim2 feature (col i is the features for the signal in column i of x)
#
# Example:
#    WaveformGen(10,0.01,1,1,0,0,0.1,'triangle');
#    WaveformGen(10,0.01,1,1,0,0,0,'square');
#
# Author: Lachlan Smith
# Work address: 8 Little Queen Street, Chippendale NSW 2008.
# email: lsmi5655@uni.sydney.edu.au
# Website: https://www.sydney.edu.au/engineering
# Janurary 2021; Last revision: 14-1-2021
#------------ BEGIN CODE --------------
def getzcfeat( x , deadzone , winsize , wininc , datawin = None ):
    
    size = np.shape(x)
    datasize = size[0]
    
    
    try:
        Nsignals = size[1]
    except:
        Nsignals = 1
    
    #check if window shape has been included, if not
    #default to rectangular window :)
    if datawin is None:
        datawin = np.ones(winsize)
        
    #print(datawin)
    
    
    numwin = math.floor(((datasize - winsize) / wininc ) + 1)
    
    #allocate memory
    feat = np.zeros( (numwin,Nsignals) );
    
    st = 0;
    en = winsize;
    
    
    for i in range(numwin):
        
        y = np.transpose(x[st:en]) * np.tile(datawin,(Nsignals))
        #print(y)
       
        
        
        y = ((1*( y > deadzone )) - (1*( y < -deadzone )))  
        
        # forces the zeros towards either the positive or negative
        # the filter is chosen so that the most recent +1 or -1 has
        # the most influence on the state of the zero.
        a = 1
        b = np.exp(-(np.arange(1,(winsize/2))))
        z = lfilter(b, a, y)
        
        z = (1*(z > 0)) - (1*(z < -0))    
        dz = np.diff(z)
        
        
        feat[i,:] = np.sum( 1*(np.abs(dz) == 2) )
        
        st = st + wininc
        en = en + wininc
    
    return feat



#------------- END OF CODE --------------

