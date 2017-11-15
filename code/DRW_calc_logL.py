import numpy as np 
import os
import datetime
from itertools import product
from scipy.optimize import minimize
import celerite
from celerite import terms
from astropy.table import Table
from astropy.table import vstack
from astropy.table import Column
import os 
import datetime

outDir = os.path.join(os.getcwd(),'logL',
                      datetime.datetime.now().strftime('%Y-%m-%d')+ '/')
if not os.path.exists(outDir): os.system('mkdir %s' % outDir)
print('We will save this figure in  %s'%outDir)

DirIn = 'DRWtestCeleriteZI/'
files = os.listdir(DirIn)

# Fitting  : each light curve is fit with various settings : 

sigma_in = 0.2
tau_in = 100
kernel = terms.RealTerm(log_a = 2 * np.log(sigma_in) , 
                        log_c = np.log(1/tau_in))


# Make a grid for logL...
# in units of tau / tau_in ,     
# sigma / sigma_in 
step  = 0.01
start = 0.4
stop = 2.5
N = int((stop - start) /step)
grid = np.linspace(start, stop, N)

sigma_grid = grid * sigma_in
tau_grid = grid * tau_in

log_a_grid = 2 * np.log(sigma_grid)
log_c_grid = np.log(1/tau_grid)


# initialize dictionaries to store the results  
results = {}
# define which priors we would like to try ... 
priors = ['flat', 'p1', 'p2', 'Jeff1', 'Jeff2']
for prior in  priors:
    results[prior] = {'sigma_fit':np.zeros(len(files), dtype=float),
                     'tau_fit':np.zeros(len(files), dtype=float)
                     }

for k in range(len(files)) :
    # read in saved light curves
    lc = Table.read(DirIn +files[k], format='ascii', 
                    names=['time', 'mag', 'err'] )    
    t,y,yerr = lc['time'], lc['mag'], lc['err']

    print('Light curve %s : %d / %d'%(files[k], k, len(files)))
    # call the model  with a chosen kernel instance 
    gp = celerite.GP(kernel, mean=np.mean(y))
    gp.compute(t, yerr)

    loglike_dic = {}
    for prior in priors:  
        # loop over the likelihood space .... 
        print('Calculating the negloglike grid for %s prior'%prior)
        loglike = np.zeros([N,N], dtype=float)
        for i in range(len(log_a_grid)):
            for j in range(len(log_c_grid)):
                params = [log_a_grid[i],log_c_grid[j]]

                if prior is 'flat' : 
                    def neg_log_like(params, y, gp):
                        gp.set_parameter_vector(params)
                        return -gp.log_likelihood(y)
                    
                if prior is 'p1' : # sigma*tau 
                    def neg_log_like(params,y,gp):
                        gp.set_parameter_vector(params)
                        log_a = params[0]
                        log_c = params[1]
                        return -gp.log_likelihood(y) - (log_a / 2.0) + log_c

                if prior is 'p2' : # sigma_hat * tau 
                     def neg_log_like(params, y, gp):
                        gp.set_parameter_vector(params)
                        log_a = params[0]
                        log_c = params[1]
                        return -gp.log_likelihood(y)  +0.5* (-np.log(2) + log_c - log_a  )

                if prior is 'Jeff1' : # (1/sigma) * (1/tau) 
                    def neg_log_like(params,y,gp):
                        gp.set_parameter_vector(params)
                        log_a = params[0]
                        log_c = params[1]
                        return -gp.log_likelihood(y) + (log_a / 2.0) - log_c


                if prior is 'Jeff2' : # (1/sigma_hat) * (1/tau) 
                    def neg_log_like(params, y, gp):
                        gp.set_parameter_vector(params)
                        log_a = params[0]
                        log_c = params[1]
                        return -gp.log_likelihood(y)  +0.5* (np.log(2) - log_c + log_a  )

                # evaluate the neg_log_like on a grid of sigma, tau ...     
                loglike[i,j] = neg_log_like(params,y,gp)

        # store all logL spaces for that prior ... 
        loglike_dic[prior] = loglike

    # save in  a file ...
    fname = files[k][:-4]+'_logL.npy'
    np.save(outDir + fname, loglike_dic)
    print('Saved logL dic as %s'%fname)




