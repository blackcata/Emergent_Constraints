#!/usr/bin/env python
# coding: utf-8
#------------------------------------------------------------------------------------#
#                                                                                    #
#         SCRIPT for basic functions to calculate emergent constraint method         #
#                                             (based on the Cox et al., 2013 paper)  #
#                                                                                    #
#                                                               BY   : KM.Noh        #
#                                                               DATE : 2023.03.15    #
#                                                                                    #
#------------------------------------------------------------------------------------#

## Modules for Calculate netCDF 
import numpy    as np
import xarray   as xr

##########################################################################################
######                 Emergent Constraint Method (Cox et al 2013)
##########################################################################################

### Calculate Gaussian PDF using mean and std
def calc_GAUSSIAN_PDF(mu,sigma,x):
    
    PDF = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-1/2*((x-mu)/sigma)**2)
    
    return PDF


### Calculate PDF of Constrained Projections
def calc_PDF_EC(tmp_x,tmp_y,x,y,PDF_x):
    dx = x[1]-x[0]

    ### Inter-model Spread in CMIP Data
    xn = tmp_x.values
    yn = tmp_y.values
    
    ### Linear Regression for inter-model Diversity
    N = len(xn)
    sigma_x  = np.sqrt(np.var(xn))
    sigma_xy = np.sqrt(np.cov(xn,yn))[0,1]

    b = (sigma_xy/sigma_x)**2
    a = -1/N*(b*xn-yn).sum()

    fn = a + b*xn
    
    ### Prediction Error
    s        = np.sqrt(1/(N-2)*((yn-fn)**2).sum())
    sigma_b  = (s/sigma_x) * np.sqrt(N)

    ### Calculate PDF(y) [posteriori] by given PDF(y|x) [Priori] and PDF(x) [OBS]
    PDF_y = np.zeros(len(y))
    for ind_y in range(len(y)):
        for ind_x in range(len(x)):
            sigma_fx       = s * np.sqrt(1 + 1/N + (x[ind_x]-xn.mean())**2/(N*sigma_x**2))
            PDF_yx         = 1/np.sqrt(2*np.pi*sigma_fx**2) \
                           * np.exp(-(y[ind_y]-(a+b*x[ind_x]))**2/(2*sigma_fx**2))            
            PDF_y[ind_y] += PDF_yx * PDF_x[ind_x] * dx 
            
    thres     = 0.341 ###
    sigma_y   = find_std_from_PDF(thres,y,PDF_y)
    mean_y    = y[PDF_y.argmax()]
    
    return PDF_y, sigma_y, mean_y


### Calculate standard-deviation of Gaussian PDF
def find_std_from_PDF(thres,y,PDF):
    ind_max = PDF.argmax()

    for ind_point in range(len(y)):
        PDF_sum = PDF[ind_point:ind_max+1].sum()/PDF.sum()
        if (PDF_sum < thres): 
            sigma = y[ind_max]-y[ind_point]
            break
        
    return sigma



### Calculate the Prior of Probability in CMIP projections
def calc_PDF_EC_PRIOR(tmp_x,tmp_y,x,y):
    xn = tmp_x.values
    yn = tmp_y.values

    N = len(xn)
    sigma_x  = np.sqrt(np.var(xn))
    sigma_xy = np.sqrt(np.cov(xn,yn))[0,1]

    b = (sigma_xy/sigma_x)**2
    a = -1/N*(b*xn-yn).sum()

    fn = a + b*xn
    
    ## prediction error
    s = np.sqrt(1/(N-2)*((yn-fn)**2).sum())

    sigma_b = (s/sigma_x) * np.sqrt(N)
    sigma_fx = s * np.sqrt(1 + 1/N + (x-xn.mean())**2/(N*sigma_x**2))

    PDF_yx = 1/np.sqrt(2*np.pi*sigma_fx**2) * np.exp(-(y-(a+b*x))**2/(2*sigma_fx**2))

    return PDF_yx,sigma_fx, a+b*x