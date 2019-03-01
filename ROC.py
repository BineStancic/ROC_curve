import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.stats import poisson, binom
import scipy


#gaussian class
class gaussian:
    
    # constructor where xhigh and xlow are substitutes for +/- inf
    def __init__(self, mean, sigma):
        self.mean = mean
        self.sigma = sigma
        self.shape = lambda x: np.exp( -((x-self.mean)**2) / (2.0 * self.sigma**2 ) ) / (self.sigma*np.sqrt(2.0*np.pi ))
        self.xhigh = mean + (6.*self.sigma) 
        self.xlow = mean - (6.*self.sigma)
    
    #normalise over the range defined above
    def normalisation(self):
        return integrate.quad(self.shape, self.xlow, self.xhigh)[0]
    
    
    # integral between tcut point and xhigh (pion)
    def integralAbove(self, tcut):
        return integrate.quad(self.shape, tcut, self.xhigh)[0]
        norm = integral/self.normalisation()
        return norm
        
    # integral between xlow and tcut point (kaon)
    def integralBelow(self, tcut):
        integral = integrate.quad(self.shape, self.xlow, tcut)[0]
        norm = integral/self.normalisation()
        return norm

#ROC class
class ROC:
    
    #constructor
    def __init__(self, mom, imax, sigma):
        self.mom = mom
        self.imax = imax
        self.sigma = sigma
    
    #define the time of flight function
    def tof(self, mass):
        return(L/c * np.sqrt(1. + (mass/self.mom)**2))
        
    #Define the function that calculates TOF, from that generates a list of tcut values
    #And uses the gaussian class to find alpha and beta values
    def alpha_beta_gaussian(self):
        tof_pion = self.tof(mass_pion)
        tof_kaon = self.tof(mass_kaon)
        tcut = np.linspace(tof_pion, tof_kaon, self.imax)
        
        kaon_gaussian = gaussian(tof_kaon, self.sigma)
        pion_gaussian = gaussian(tof_pion, self.sigma)
        
        alpha = []
        beta = []
        for i in (tcut):
            alpha += [kaon_gaussian.integralBelow(i)]
            beta += [pion_gaussian.integralAbove(i)]
        return alpha, beta
		
#define parameters
L = 20.
c = 3. * 10**8
imax = 100
sigma = 4. *10**(-10)

#masses and momenta in GeV/c^2 and GeV/c respectively
mass_pion = 0.1396
mass_kaon = 0.4937
mom_1 = 3.
mom_2 = 4.
mom_3 = 6.

#unpack classes and plotting inside main
def main():  
    
    #import the ROC class and get the alpha and beta arguments for all 3 mom
    mom1 = ROC(mom_1, imax, sigma)
    alpha, beta = mom1.alpha_beta_gaussian()
    plt.plot(alpha, beta, label = 'mom = 3GeV/c')

    mom2 = ROC(mom_2, imax, sigma)
    alpha, beta = mom2.alpha_beta_gaussian()
    plt.plot(alpha, beta, label = 'mom = 4GeV/c')

    mom3 = ROC(mom_3, imax, sigma)
    alpha, beta = mom3.alpha_beta_gaussian()
    plt.plot(alpha, beta, label = 'mom = 6GeV/c')

    plt.xlabel('alpha')
    plt.ylabel('beta')
    plt.legend()
    plt.show()

main()