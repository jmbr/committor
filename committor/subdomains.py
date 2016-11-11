__all__ = ['make_subdomains']

import ast
import configparser
from typing import Tuple

import numpy as np

import fenics


class InvalidConfig(Exception):
    """Invalid configuration file."""
    pass


def make_subdomains(cfg_file: str) \
        -> Tuple[fenics.CompiledSubDomain, fenics.CompiledSubDomain, int]:
    """Create reactant and product subdomains.

    This function reads a specification of the reactant (R) and the
    product (P) from a configuration file and generates two
    (optimized) instances of the SubDomain class to identify R and P.

    Parameters
    ----------
    cfg_file : str
        Configuration file containing the descriptions of the reactant
        and the product. The file must contain two sections named
        'reactant' and 'product' defines as d-dimensional spheres by
        specifying their center and radius.

    Returns
    -------
    reactant : fenics.CompiledSubDomain
        Instances of CompiledSubDomain suitable for use in the
        committor equation solver.

    product : fenics.CompiledSubDomain
        Instances of CompiledSubDomain suitable for use in the
        committor equation solver.

    dim : int
        Dimension of the subdomains.

    Raises
    ------
    InvalidConfig
        In case the configuration file is ill-formed.

    """
    parser = configparser.ConfigParser()
    parser.read(cfg_file)

    sections = {section.lower() for section in parser.sections()}
    if 'reactant' not in sections:
        raise InvalidConfig('Unable to find definition of reactant in {!r}'
                            .format(cfg_file))
    if 'product' not in sections:
        raise InvalidConfig('Unable to find definition of product in {!r}'
                            .format(cfg_file))

    dim = None
    subdomains = {}
    for domain in parser.sections():
        section = parser[domain]

        center, radius = parse_center_and_radius(section, dim)
        dim = center.shape[0]

        subdomains[domain.lower()] = generate_subdomain(center, radius)

    return subdomains['reactant'], subdomains['product'], dim


def parse_center_and_radius(section: configparser.SectionProxy, dim: int = None) \
        -> Tuple[np.array, float]:
    """Parse center and radius from a configuration file section.

    """
    center = np.array(ast.literal_eval(section['center']), dtype=np.float)
    radius = float(section['radius'])

    if dim is not None and center.shape[0] != dim:
            raise InvalidConfig('Inconsistent number of dimensions')

    if radius <= 0.0:
        raise InvalidConfig('Radius is non-positive')

    return center, radius


def generate_subdomain(center: np.array, radius: float) \
        -> fenics.CompiledSubDomain:
    """Generate spherical subdomain.

    """
    code = []
    for i in range(center.shape[0]):
        code.append('(x[{0}] - ({1})) * (x[{0}] - ({1}))'
                    .format(i, center[i]))
    code = ' + '.join(code)

    cppcode = code + ' < {} + DOLFIN_EPS'.format(radius**2)

    return fenics.CompiledSubDomain(cppcode)
