import numpy as np

class _ModelGeneral(object):
    """
    Code that is common to different Model classes
    """
    def _init_finish(self):
        """common parts of __init__()"""
        self.n_params = len(self._params)

    def set_values(self, values):
        """
        Set the values of parameters

        Input :
            Values: *np.array*
                Vector of values to be set
        """
        if len(values) != self.n_params:
            msg = 'wrong number of params: {:} vs. {:}'
            raise ValueError(msg.format(len(values), self.n_params))

        for (parameter, value) in zip(self._params_names, values):
            self.params[parameter] = value


class _ModelGeneralBrokenPowerLaw(_ModelGeneral):
    def _get_log_prior(self):
        """check if A and q_br are positive"""
        if self._params['A'] > 0 and self._params['q_br'] > 0:
            return 0.0
        else:
            return -np.inf

class Model4Parameter(Model5Parameter):
    """
    f() that has 1 slope in s and 2 slopes in q with q_br fixed at 1.7e-4
    """
    def __init__(self):
        self._params_names = ['A', 'm', 'n', 'p']
        self._params_latex = ['$A$', '$m$', '$n$', '$p$']

        self._params['q_br'] = 1.7e-4


class Model5Parameter(_ModelGeneralBrokenPowerLaw):
    """
    f() that has 1 slope in s and 2 slopes in q
    """
    def __init__(self):
        self._params_names = ['A', 'm', 'n', 'p', 'q_br']
        self._params_latex = ['$A$', '$m$', '$n$', '$p$', '$q_\mathrm{br}$']


