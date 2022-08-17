"""Functions used in the simulation which cannot be represented by simple analytical functions.

These functions are handled here because they cannot be written in a form that DAE Tools knows how
to automatically differentiate. For example, they may contain `if` statements or a function from an
external library that the DAE Tools library doesn't know about.
"""
import scipy.interpolate as sintrp

from daetools.pyDAE import daeScalarExternalFunction, adouble


class InterpTimeScalar(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, time, tvec, yvec):
        arguments = {}
        arguments["time"] = time
        self.cache = None
        self.interp = sintrp.interp1d(tvec, yvec, bounds_error=False,
                                      fill_value=yvec[-1])
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
