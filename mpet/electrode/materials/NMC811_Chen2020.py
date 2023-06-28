import numpy as np


def NMC811_Chen2020(self, y, ybar, T, muR_ref):
    """
    This function was obtained from Chen et al. 2020.
    """
    OCV = (-0.8090*y + 4.4875 - 0.0428*np.tanh(18.5138*(y - 0.5542)) 
           - 17.7326*np.tanh(15.7890*(y - 0.3117)) 
           + 17.5842*np.tanh(15.9308*(y - 0.3120)))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR