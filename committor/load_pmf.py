#!/usr/bin/env python

import sys
from typing import Tuple

import numpy as np
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt

import fenics


class DimensionsExceeded(Exception):
    """Attempt to use a space of dimension higher than three."""
    def __init__(self, dim: int) -> None:
        self.dim = dim

    def __str__(self):
        return 'Maximum number of dimensions exceeded: {} > 3'.format(self.dim)


def get_line_or_fail(f):
    """Obtain a line and verify that it is well formed."""
    line = f.readline()
    if not line.startswith('#'):
        print('Error: Invalid PMF file {!r}.'.format(file_name),
              file=sys.stderr)
        sys.exit(-1)
    
    return line
    

def get_num_variables(f):
    """Parse number of variables from open file."""
    line = get_line_or_fail(f)

    _, num_variables_str = line.split()

    return int(num_variables_str)


def get_variable_data(f):
    line = get_line_or_fail(f)

    _, minimum, width, size, pbc = line.split()

    return int(minimum), int(width), int(size), pbc == 1



def load_pmf_file(pmf_file_name: str) -> Tuple[np.array, np.array]:
    """Load potential of mean force from file."""
    data = np.loadtxt(pmf_file_name)

    x = data[:, 0:-1]
    F = data[:, -1]

    return x, F


def load_est_file(pmf_file_name: str) -> Tuple[np.array, np.array]:
    """Load estimated gradient of the potential of mean force from file."""
    data = np.loadtxt(pmf_file_name)

    n = int(data.shape[1] / 2)
    x = data[:, 0:n]
    dF = data[:, n:]

    return x, dF


def load_pmf(grad_file_name, pmf_file_name, est_file_name):
    with open(grad_file_name) as f:
        num_variables = get_num_variables(f)
        if num_variables > 3:
            raise DimensionsExceeded(num_variables)

        scalar_dim = 1

        minima = []
        maxima = []
        widths = []
        sizes = []
        pbcs = []

        for n in range(num_variables):
            minimum, width, size, pbc = get_variable_data(f)
            maximum = minimum + width * size

            minima.append(minimum)
            maxima.append(maximum)
            widths.append(width)
            sizes.append(size)
            pbcs.append(pbc)

            scalar_dim *= size

    x, F = load_pmf_file(pmf_file_name)
    assert x.shape[0] == scalar_dim

    x, dF = load_est_file(est_file_name)
    assert x.shape[0] == scalar_dim

    return x, F, dF
        


def main(args):
    if len(args) < 4:
        print('Usage:', sys.argv[0], 'GRAD-FILE PMF-FILE EST-FILE',
              file=sys.stderr)
        sys.exit(-1)

    x, F, dF = load_pmf(args[1], args[2], args[3])

    xx, yy = np.mgrid[-180:180, -180:180]
    dFdx = scipy.interpolate.NearestNDInterpolator(x, dF[:, 0])
    zz = dFdx(xx, yy)
    # zz = scipy.interpolate.griddata(x, dF[:, 0], (xx, yy), method='linear')



    ax1 = plt.subplot(121)
    # ax1.pcolormesh(xx, yy, zz, cmap='RdBu_r')
    ax1.imshow(zz.T, extent=(-180, 180, -180, 180), origin='lower',
               cmap='RdBu_r')
    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-180, 180])
    ax1.set_xticks(np.linspace(-180, 180, 5))
    ax1.set_yticks(np.linspace(-180, 180, 5))
    ax1.set_aspect('equal')

    ax2 = plt.subplot(122)
    n = int(np.sqrt(x.shape[0]))
    ax2.pcolormesh(dF[:, 0].reshape((n, n)).T, cmap='RdBu_r')
    # ax2.imshow(dF[:, 0].reshape((n, n)).T, extent=(-180, 180, -180, 180),
    #            origin='lower', cmap='RdBu_r')
    ax2.set_xlim([0, 72])
    ax2.set_ylim([0, 72])
    # ax2.set_xticks(np.linspace(-180, 180, 5))
    # ax2.set_yticks(np.linspace(-180, 180, 5))
    ax2.set_aspect('equal')
    plt.show()

    # zzz = scipy.interpolate.bisplev(xx[:, 0], yy[0, :], tck)


    # # fenics.set_log_level(fenics.DBG)
    # 
    # if num_variables == 1:
    #     mesh = fenics.IntervalMesh(sizes[0], minima[0], maxima[0])
    # elif num_variables == 2:
    #     p0 = fenics.Point(*minima)
    #     p1 = fenics.Point(*maxima)
    #     mesh = fenics.RectangleMesh(p0, p1, *sizes)
    # elif num_variables == 3:
    #     p0 = fenics.Point(*minima)
    #     p1 = fenics.Point(*maxima)
    #     mesh = fenics.BoxMesh(p0, p1, *sizes)
    # 
    # V = fenics.FunctionSpace(mesh, 'Lagrange', 2)
    # 
    # bc = fenics.DirichletBC(V, fenics.Constant(0), boundary)
    # 
    # u = fenics.TrialFunction(V)
    # v = fenics.TestFunction(V)
    # 
    # a = fenics.dot(fenics.grad(u), fenics.grad(v)) * fenics.dx
    # L = fenics.Constant(0) * fenics.dx
    # 
    # u = fenics.Function(V)
    # fenics.solve(a == L, sol, bc)

    
if __name__ == '__main__':
    main(sys.argv)


# main(['', 'ala3.grad', 'ala3.grad.pmf', 'ala3.grad.est'])


# print(mesh.coordinates())

# fenics.plot(mesh)

# data = np.loadtxt(est_file_name)
# x = data[:, 0:num_variables]
# dF = data[:, num_variables:]

# plot(entries, scalar_dim, num_variables, x, F, dF)
