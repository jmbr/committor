import string

import fenics

from potential_mean_force import PotentialMeanForce


__all__ = ['generate_periodic_boundary_conditions']


CPPCODE = '''
class PeriodicBoundaryConditions : public SubDomain {
 public:
  bool inside(const Array<double>& x, bool on_boundary) const {
    return ((${inside}) && on_boundary);
  }

  void map(const Array<double>& x, Array<double>& y) const {
    for (std::size_t i = 0; i < x.size(); ++i)
      y[i] = x[i];

${maps}
  }
};
'''


def generate_periodic_boundary_conditions(PMF: PotentialMeanForce) \
        -> fenics.CompiledSubDomain:
    """Generate periodic boundary conditions on the fly.

    Automatically generate C++ code to implement periodic boundary
    conditions on the relevant collective variables.

    Parameters
    ----------
    PMF : PotentialMeanForce
        Potential of mean force object containing all the information
        about the space of collective variables.

    Returns
    -------
    subdomain : fenics.CompiledSubDomain
        Compiled code to identify and wrap the periodic boundaries.

    """
    minima = PMF.minima
    maxima = PMF.maxima

    inside_expr = []
    maps_expr = []

    for n in range(PMF.num_variables):
        if not PMF.pbcs[n]:
            continue

        inside_expr.append('std::fabs(x[{}] - ({})) < DOLFIN_EPS'
                           .format(n, minima[n]))

        maps_expr.append('    if (std::fabs(x[{0}] - ({1})) < DOLFIN_EPS)\n'
                         '      y[{0}] -= {2};'
                         .format(n, maxima[n], maxima[n] - minima[n]))

    inside_str = ' || '.join(inside_expr)
    maps_str = '\n'.join(maps_expr)

    template = string.Template(CPPCODE)
    code = template.safe_substitute(inside=inside_str, maps=maps_str)

    return fenics.CompiledSubDomain(code)
