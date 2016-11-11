import fenics

from potential_mean_force import PotentialMeanForce
from periodic_boundary_conditions import generate_periodic_boundary_conditions
import subdomains


class CommittorError(Exception):
    """Error computing the committor function."""
    pass


def generate_mesh(PMF: PotentialMeanForce, num_cells: int) -> fenics.Mesh:
    """Generate a mesh for the space of collective variables.

    """
    p0 = fenics.Point(*PMF.minima)
    p1 = fenics.Point(*PMF.maxima)

    n = PMF.num_variables
    if n == 1:
        mesh = fenics.IntervalMesh(num_cells, PMF.minima[0], PMF.maxima[0])
    elif n == 2:
        mesh = fenics.RectangleMesh(p0, p1, num_cells, num_cells)
    elif n == 3:
        mesh = fenics.BoxMesh(p0, p1, num_cells, num_cells, num_cells)
    else:
        raise CommittorError('Invalid number of dimensions: {}'.format(n))

    return mesh


def committor(grad_file: str, pmf_file: str, est_file: str, output: str,
              num_cells: int, cfg_file: str, verbose: bool) -> None:
    """Compute committor function.

    Parameters
    ----------
    grad_file : str
        File name of gradient from colvars.

    pmf_file : str
        File name containing the potential of mean force file from
        colvars.

    est_file : str
        File name of the estimated gradient of the potential of mean
        force from colvars.

    output_file : str
        File name where to save results.

    num_cells : int
        Number of cells per dimension in the mesh. This controls the
        accuracy of the solution.

    cfg_file : str
        Configuration file where the reactant and the product are
        defined.

    verbose : bool
        Whether to print verbose progress information or not.

    Raises
    ------
    CommittorError
        Error computing the committor function.

    """
    fenics.set_log_level(fenics.PROGRESS if verbose else fenics.INFO)

    PMF = PotentialMeanForce(grad_file, pmf_file, est_file)
    PMF.load()

    mesh = generate_mesh(PMF, num_cells)

    PBC = generate_periodic_boundary_conditions(PMF)

    V = fenics.FunctionSpace(mesh, 'Lagrange', 1, constrained_domain=PBC)

    u = fenics.TrialFunction(V)
    v = fenics.TestFunction(V)

    reactant, product, subdomain_dim = subdomains.make_subdomains(cfg_file)
    if subdomain_dim != PMF.num_variables:
        raise CommittorError('Inconsistent subdomain dimension.')

    bc_reactant = fenics.DirichletBC(V, fenics.Constant(0.0), reactant)
    bc_product = fenics.DirichletBC(V, fenics.Constant(1.0), product)
    dirichlet_boundaries = [bc_reactant, bc_product]

    dF = PMF.get_gradient(degree=1)

    a = (fenics.dot(fenics.grad(u), fenics.grad(v)) * fenics.dx
         - fenics.dot(fenics.grad(u), dF) * v * fenics.dx)
    L = fenics.Constant(0) * v * fenics.dx

    u = fenics.Function(V)

    fenics.solve(a == L, u, dirichlet_boundaries)

    vtkfile = fenics.File(output)
    vtkfile << u
