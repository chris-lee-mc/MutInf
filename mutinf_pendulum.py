#!~/usr/bin/python

PI = 3.141592653589793238462643383279502884197

from numpy import *
from scipy import integrate as integrate

r1 = 5.885
r2 = 5.885
q1 = 1.0
q2 = 1.0
vacuum_permittivity = 8.8541878160E-12
echarge = 1.602177330E-19
boltzmann = 1.3806580E-23
temperature = 298.0
avagadro = 6.02213670E+23
kT_in_joules = boltzmann * temperature *avagadro
print kT_in_joules / 4.184
echarge2_over_vacuum = echarge * echarge / (vacuum_permittivity * 4 * PI)
print "energy for division by angstrom distance, in Ang:"+str(((avagadro*echarge2_over_vacuum/(1.0 * kT_in_joules)) / (1E-10 )))
class pendulum:
    r = 5.885E-10
    q = 1.0
    eps = 4.0
    def __init__(self, r, q):
        self.r = r * 1.0E-10
        self.q = q
#    def __init__(r, q, eps):
#        self.r = r * 1.0E-10
#        self.q = q
#        self.eps = eps
#    def __init__():
#        self.r = 5.885E-10
#        self.q = 1.0
#    def __init__(q):
#        self.r = 5.885E-10 
#        self.q = q


class twobody(pendulum):
    partner = []
    L = []
    eps = 4.0
    Z = None
    def __init__(self, r1, q1, r2, q2, eps1, eps2, L):
        self.r = r1 * 1.0E-10 #Angstroms to meters
        self.q = q1 
        self.partner = pendulum(r2, q2)
        self.L = L * 1.0E-10 #Angstroms to meters
        self.eps = max(eps1, eps2)
        self.Z = self.__partition_func__()
        print "test energy at 0,pi ="+str(avagadro*self.energy(0,PI)/kT_in_joules)
    def myzero(self,x):
        return 0.0
    def mypitwothirds(self,x):
        return 2.0*PI/3.0
    def my2pi(self,x):
        return 2.0*PI
    def mypi(self,x):
        return PI
    def mynpi(self,x):
        return -PI
    def nmypiover2(self,x):
        return -PI/2.0


    def energy(self, mtheta, ptheta):
        return (echarge2_over_vacuum * self.q * self.partner.q / (self.eps * sqrt((self.L + self.partner.r * cos(ptheta) - self.r * cos(mtheta))**2 + (self.partner.r * sin(ptheta) - self.r * sin(mtheta))**2)))
        

    def boltz_factor(self, ptheta, mtheta): #note arguments reversed for dblquad
        return (exp(avagadro*-self.energy(mtheta,ptheta)/kT_in_joules))

    def __partition_func__(self):
        partfunc = integrate.dblquad(self.boltz_factor, 0, PI/2.0, self.nmypiover2, self.myzero, epsabs=1.0E-12, epsrel=1.0E-12)
        print "partition function = "+str(partfunc)
        return partfunc[0]

    def MI_term(self, ptheta, mtheta):
        return (1.0/self.Z) * self.boltz_factor(ptheta,mtheta) * log(0.25*PI*(PI)*(1.0/self.Z)*self.boltz_factor(ptheta,mtheta))

    def mutual_info(self):
        return integrate.dblquad(self.MI_term, 0.0, PI/2.0, self.nmypiover2, self.myzero, epsabs=1.0E-12, epsrel=1.0E-12)

for x in range(20):
    test = twobody(5.885,1.0,5.885,1.0, 4.0, 4.0, 6+x)
    #print test.MI_term(0,0)
    print 6+x, ":",test.mutual_info()[0]


#END
        
    
