from pyomo.environ import *
import numpy as np

_Mean_Penalty=0
_Var_Penalty=1.0

def mk_objective(model, kappa, Z, V):
    
    def getA_term(model, kappa, Z, V, m, i):
        return model.b[m] * 2 * Z[i] ** 2 * (Z[m + 1] - Z[m]) / (Z[m + 1] * Z[m]) / kappa[m]

    def getB_term(kappa, Z, V, m):
        return (Z[m + 1] - Z[m]) / (Z[m + 1] * Z[m])
    

    terms = []
    for i in range(1, kappa.size + 1):
        B = (V[0] - sum(getB_term(kappa, Z, V, m) for m in range(0, i))) * (Z[i] ** 2)
        terms += [2 * getA_term(model, kappa, Z, V, m, i) * (B - V[i]) for m in range(0, i) ]              
        terms += [getA_term(model, kappa, Z, V, m, i) * getA_term(model, kappa, Z, V, n, i) for m in range(0, i) for n in range(0, i)]
    return sum(terms)

def compute_rates(data, max_step_size):
    """
    Compute bifurcation and annihilation rates from summary statistics.

    The estimator expects Sholl-plot summary statistics (mean and variance per
    radial bin) and summary statistics of bifurcation counts (mean and variance).
    These quantities are used to infer the event rates of a branching-and-
    annihilating process.

    Parameters
    ----------
    data : dict
      Input container with the following structure:

      data = {
        "sholl": {
          "bin_size": float,
          "mean": numpy.ndarray,   # shape (K,)
          "var":  numpy.ndarray,   # shape (K,)
        },
        "bifurcations": {
          "mean": float,
          "var":  float,
        },
      }

      Where:
      - `data["sholl"]["bin_size"]` is the spatial bin size used to build the Sholl plot
      - `data["sholl"]["mean"][i]` is the mean Sholl intersection count in bin i
      - `data["sholl"]["var"][i]` is the variance of the Sholl intersection count in bin i
      - `data["bifurcations"]["mean"]` is the mean number of bifurcations
      - `data["bifurcations"]["var"]` is the variance of the number of bifurcations

    max_step_size : float
      Maximum advancement (in distance from the soma) allowed for a single
      elongation step in the model. This value bounds the radial increment used
      by the estimator and should be expressed in the same spatial units as the
      Sholl binning.

    Returns
    -------
    dict
      Dictionary containing the estimated rates and any additional derived values
      produced by the implementation. At minimum, the returned dictionary is
      expected to include:

      - "bifurcation_rate"
      - "annihilation_rate"

    Notes
    -----
    - `data["sholl"]["mean"]` and `data["sholl"]["var"]` must be 1D arrays of equal length
    - Variances must be non-negative
    - Ensure `bin_size` and `max_step_size` use consistent spatial units
    """    

##    Solves a QP problem using Pyomo with vector-style variable indexing.
##
##    Minimize: 0.5*x[0]^2 + x[1]^2 + x[0]*x[1] + 3*x[0]
##    Subject to: x[0] + x[1] >= 1, x[i] >= 0
##
##    Returns:
##        np.ndarray: [x[0], x[1], objective_value]


    global _Mean_Penalty, _Var_Penalty
    dx = data['sholl']['bin_size']
    Z =  data['sholl']['mean']
    V = data['sholl']['var']

    if 'bifurcations' in data:
        n_bif = [data['bifurcations']['mean'], data['bifurcations']['var']]
    else:
        n_bif = None

    # get the kappa
    kappa = np.log(Z[1:] / Z[:-1]) / dx

    if np.any(kappa == 0):
        kappa += 1e-5
        
    model = ConcreteModel()

    # Define index set and variables
    model.b = Var(range(kappa.size), domain=NonNegativeReals)

    # define 1 slack variables for eventual constraints of variance of bifurcations
    model.s = Var(range(2), domain=Reals)
            
    # Constraint: 
    model.constraints = ConstraintList()
    for i in range(kappa.size):
        model.constraints.add(model.b[i] >= kappa[i])
        model.constraints.add((2 * model.b[i] - kappa[i]) * max_step_size <= 1)

    # if we have number of bifurcations as contraings
    if n_bif:        
        # constraint the average number of bifurcations
        if n_bif[0]:
            f = (Z[1:] - Z[:-1]) / kappa
            model.constraints.add(sum(f[i] * model.b[i] for i in range(kappa.size)) == n_bif[0])
        
        # constrain the variance for the number of bifurcations
        if n_bif[1]:            
            f1 = (Z[1:] - Z[:-1]) / kappa
            f2 = - Z[:-1] * (np.power(Z[1:] / Z[:-1], 2) - 2 * kappa * dx * Z[1:] / Z[:-1] - 1) / kappa ** 2 + V[:-1] * np.power((Z[1:] - Z[:-1]) / Z[:-1] / kappa, 2)
            f3 = 2 * Z[:-1] * (np.power(Z[1:] / Z[:-1], 2) - 2 * kappa * dx * Z[1:] / Z[:-1] - 1) / kappa ** 3
            
            var_terms = [ ]
            for i in range(kappa.size):
                var_terms.append(f1[i] * model.b[i] + f2[i] * model.b[i] ** 2 + f3[i] * model.b[i] ** 3)

            for i in range(1, kappa.size):
                for j in range(i + 1, kappa.size + 1):
                    term1 = model.b[i - 1] * (Z[i] - Z[i - 1])/  Z[i - 1] / kappa[i - 1] 
                    term2 = model.b[j - 1] * (Z[j] - Z[j - 1]) / Z[j - 1] / kappa[j - 1]
                    for k in range(1, i + 1):
                        if k == 1:
                            Vprev = V[0]
                        term3 = (2 * model.b[k - 1] - kappa[k - 1]) / kappa[k - 1] * Z[k] * (Z[k] - Z[k - 1]) / Z[k - 1] + Vprev * np.power(Z[k] / Z[k - 1], 2)
                        var_terms.append(2 * term1 * term2 * Z[j - 1] / Z[i - 1] * term3)
                        Vprev = term3
                
            model.constraints.add(sum(var_terms) + model.s[1] == n_bif[1]) 

    # Objective
    model.obj = Objective(
        expr=mk_objective(model, kappa, Z, V) + _Mean_Penalty * model.s[0] ** 2 + _Var_Penalty * model.s[1] ** 2,
        sense=minimize
    )

    
    # Solve
    solver = SolverFactory('ipopt')
    #solver.options['max_iter'] = max_iter
    solver.solve(model, tee=False)
    #print("Objective value:", value(model.obj))
    #print(value(model.s))
    
    # Extract bifurcation rates as array
    b = np.array([value(model.b[i]) for i in model.b])
    b[b < 0] = 0.

    # calculate annihilation rates
    a = - kappa + b
    a[a < 0] = 0.
    return { 'bifurcation_rate':b, 'annihilation_rate':a }
