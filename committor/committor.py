import math
import sys

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

import fenics

import load_pmf


EPS = sys.float_info.epsilon


class Reactant(fenics.SubDomain):
    center = np.array([160, 120])
    radius = 20

    def inside(self, x, on_boundary) -> bool:
        z = x - self.center
        return np.dot(z, z) - self.radius**2 < 0.0


class Product(fenics.SubDomain):
    center = np.array([-70, -70])
    radius = 20

    def inside(self, x, on_boundary) -> bool:
        z = x - self.center
        return np.dot(z, z) - self.radius**2 < 0.0


class PeriodicBoundary(fenics.SubDomain):
    def inside(self, x, on_boundary):
        return ((math.isclose(x[0], -180.0, rel_tol=EPS)
                 or math.isclose(x[1], -180.0, rel_tol=EPS))
                and on_boundary)

    def map(self, x, y):
        if math.isclose(x[0], 180.0, rel_tol=EPS):
            y[0] = x[0] - 360.0
            y[1] = x[1]
        else:
            y[0] = x[0]
            y[1] = x[1] - 360.0


class PotentialMeanForce(fenics.Expression):
    dFdx = None
    dFdy = None

    def load_pmf(self):
        x, pmf, dpmf = load_pmf.load_pmf('ala3.grad',
                                         'ala3.grad.pmf',
                                         'ala3.grad.est')
        self.x = x
        self.pmf = pmf
        self.dpmf = dpmf

        self.dFdx = scipy.interpolate.NearestNDInterpolator(x, dpmf[:, 0])
        self.dFdy = scipy.interpolate.NearestNDInterpolator(x, dpmf[:, 1])

    def eval(self, values, x):
        if self.dFdx is None or self.dFdy is None:
            self.load_pmf()

        values[0] = -self.dFdx(x)
        values[1] = -self.dFdy(x)

    def value_shape(self):
        return 2,


def main():
    fenics.set_log_level(fenics.PROGRESS)

    n = 256
    
    p0 = fenics.Point(-180, -180)
    p1 = fenics.Point(180, 180)
    mesh = fenics.RectangleMesh(p0, p1, n, n, 'left/right')
    # mesh = fenics.UnitSquareMesh(n, n)
    V = fenics.FunctionSpace(mesh, 'Lagrange', 1,
                             constrained_domain=PeriodicBoundary())

    u = fenics.TrialFunction(V)
    v = fenics.TestFunction(V)

    reactant = Reactant()
    bcR = fenics.DirichletBC(V, fenics.Constant(0), reactant)
    
    product = Product()
    bcP = fenics.DirichletBC(V, fenics.Constant(1), product)

    dF = PotentialMeanForce(degree=1)

    a = (fenics.dot(fenics.grad(u), fenics.grad(v)) * fenics.dx
         - fenics.dot(fenics.grad(u), dF) * v * fenics.dx)
    L = fenics.Constant(0) * v * fenics.dx
    
    u = fenics.Function(V)
    
    fenics.solve(a == L, u, [bcR, bcP])
    
    # problem = fenics.LinearVariationalProblem(a, L, u, [bcR, bcP])
    # M = u * fenics.dx
    # tol = 1e-3
    # solver = fenics.AdaptiveLinearVariationalSolver(problem, M)
    # solver.solve(tol)
    # solver.summary()
    
    vtkfile = fenics.File('committor.pvd')
    vtkfile << u
    
    fenics.plot(u, cmap='RdBu_r')
    plt.show()


if __name__ == '__main__':
    main()


# main()
