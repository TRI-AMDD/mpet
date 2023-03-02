import numpy as np


def LiFePO4_testfit(self, y, ybar, T, muR_ref):
    muRtheta = -self.eokT*3.422
    coeff = np.array([-0.08237065703940517, 3.6048303082624367, 
                     -44.62877035306038, 297.23159214333754, 
                     -1128.1280737601173, 2371.5640960425603, 
                     -2387.1787175524187, 88.70290517794763,
                     2165.3467338166643, -1895.2144080547582, 
                     528.8593491458244])  
    muRhomog = (coeff[0]*y**0 + coeff[1]*y**1 
                + coeff[2]*y**2 + coeff[3]*y**3 
                + coeff[4]*y**4 + coeff[5]*y**5 
                + coeff[6]*y**6 + coeff[7]*y**7 
                + coeff[8]*y**8 + coeff[9]*y**9 
                + coeff[10]*y**10)
    
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
