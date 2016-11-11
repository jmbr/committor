from typing import Tuple, io

import numpy as np
import scipy.interpolate

import fenics


__all__ = ['DimensionsExceeded', 'InvalidPMFFile', 'PotentialMeanForce']


class DimensionsExceeded(Exception):
    """Attempt to use a space of dimension higher than three.

    """
    def __init__(self, dim: int) -> None:
        self.dim = dim

    def __str__(self):
        return ('Maximum number of dimensions exceeded: {} > 3'
                .format(self.dim))


class InvalidPMFFile(Exception):
    """Malformed PMF file."""
    pass


def get_line(f: io) -> str:
    """Obtain a line and verify that it is well formed.

    """
    line = f.readline()
    if not line.startswith('#'):
        raise InvalidPMFFile()

    return line


def get_num_variables(f: io) -> int:
    """Parse number of variables from open file.

    """
    line = get_line(f)

    _, num_variables_str = line.split()

    return int(num_variables_str)


def get_variable_data(f: io) -> Tuple[int, int, int, bool]:
    """Parse collective variable information.

    """
    line = get_line(f)

    _, minimum, width, size, pbc = line.split()

    return float(minimum), float(width), int(size), int(pbc) == 1


def load_pmf_file(pmf_file_name: str) -> Tuple[np.array, np.array]:
    """Load potential of mean force from file.

    """
    data = np.loadtxt(pmf_file_name)

    x = data[:, 0:-1]
    F = data[:, -1]

    return x, F


def load_est_file(pmf_file_name: str) -> Tuple[np.array, np.array]:
    """Load estimated gradient of the potential of mean force from file.

    """
    data = np.loadtxt(pmf_file_name)

    n = int(data.shape[1] / 2)
    x = data[:, 0:n]
    dF = data[:, n:]

    return x, dF


class PotentialMeanForce:
    """Potential of mean force.

    """
    def __init__(self, grad_file: str, pmf_file: str, est_file: str) \
            -> None:
        self.grad_file = grad_file
        self.pmf_file = pmf_file
        self.est_file = est_file
        self.scalar_dim = None
        self.minima = []
        self.maxima = []
        self.widths = []
        self.sizes = []
        self.pbcs = []
        self.num_variables = None
        self.x = None
        self.F = None
        self.dF = None

    def load(self) -> None:
        """Load the potential of mean force and its gradient.

        """
        with open(self.grad_file) as f:
            self.num_variables = get_num_variables(f)
            if self.num_variables > 3:
                raise DimensionsExceeded(self.num_variables)

            self.scalar_dim = 1

            for n in range(self.num_variables):
                minimum, width, size, pbc = get_variable_data(f)
                maximum = minimum + width * size

                self.minima.append(minimum)
                self.maxima.append(maximum)
                self.widths.append(width)
                self.sizes.append(size)
                self.pbcs.append(pbc)

                self.scalar_dim *= size

        self.x, self.F = load_pmf_file(self.pmf_file)
        assert self.x.shape[0] == self.scalar_dim

        self.x, self.dF = load_est_file(self.est_file)
        assert self.x.shape[0] == self.scalar_dim

    def get_gradient(self, *args, **kwargs):
        n = self.num_variables
        if n == 1:
            return GradPotentialMeanForce1D(self, *args, **kwargs)
        elif n == 2:
            return GradPotentialMeanForce2D(self, *args, **kwargs)
        elif n == 3:
            return GradPotentialMeanForce3D(self, *args, **kwargs)


class GradPotentialMeanForce1D(fenics.Expression):
    __dFdx = None

    def __init__(self, pmf: PotentialMeanForce, *args, **kwargs) -> None:
        x, d = pmf.x, pmf.dF
        self.__dFdx = scipy.interpolate.NearestNDInterpolator(x, d[:, 0])

    def eval(self, values, x):
        values[0] = -self.__dFdx(x)

    def value_shape(self):
        return 1,


class GradPotentialMeanForce2D(fenics.Expression):
    __dFdx = None
    __dFdy = None

    def __init__(self, pmf: PotentialMeanForce, *args, **kwargs) -> None:
        x, d = pmf.x, pmf.dF
        self.__dFdx = scipy.interpolate.NearestNDInterpolator(x, d[:, 0])
        self.__dFdy = scipy.interpolate.NearestNDInterpolator(x, d[:, 1])

    def eval(self, values, x):
        values[0] = -self.__dFdx(x)
        values[1] = -self.__dFdy(x)

    def value_shape(self):
        return 2,


class GradPotentialMeanForce3D(fenics.Expression):
    __dFdx = None
    __dFdy = None
    __dFdz = None

    def __init__(self, pmf: PotentialMeanForce, *args, **kwargs) -> None:
        x, d = pmf.x, pmf.dF
        self.__dFdx = scipy.interpolate.NearestNDInterpolator(x, d[:, 0])
        self.__dFdy = scipy.interpolate.NearestNDInterpolator(x, d[:, 1])
        self.__dFdz = scipy.interpolate.NearestNDInterpolator(x, d[:, 2])

    def eval(self, values, x):
        values[0] = -self.__dFdx(x)
        values[1] = -self.__dFdy(x)
        values[2] = -self.__dFdz(x)

    def value_shape(self):
        return 3,
