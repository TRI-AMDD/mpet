import numpy as np
import scipy.interpolate as sintrp

from daetools.pyDAE import daeScalarExternalFunction, adouble

class InterpScalar(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, time, xvec, yvec):
        arguments = {}
        arguments["time"] = time
        self.cache = None
        self.interp = sintrp.interp1d(xvec, yvec)
        daeScalarExternalFunction.__init__(self, Name, Model, units, arguments)

    def Calculate(self, values):
        time = values["time"]
        # A derivative for Jacobian is requested - always return 0.0
        if time.Derivative != 0:
            return adouble(0)
        # Store the previous time value to prevent excessive
        # interpolation.
        if self.cache:
            if self.cache[0] == time.Value:
                return adouble(float(self.cache[1]))
        yval = float(self.interp(time.Value))
        self.cache = (time.Value, yval)
        return adouble(yval)
