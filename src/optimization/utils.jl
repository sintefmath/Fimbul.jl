function calibrate_case(objective, case, n_steps, opt_config; lbfgs_args = NamedTuple())

    case = deepcopy(case)
    _, cfg = setup_reservoir_simulator(case;
        relaxation=true,
        tol_cnv=1e-6,
        tol_cnv_well=1e-5,
        info_level=-1,
    );
    push!(cfg[:timestep_selectors], TimestepSelector(min = convert_to_si(0.5, :day)))
    cfg[:tolerances][:Reservoir][:default] = 1e-6
    for well in well_symbols(case.model)
        cfg[:tolerances][well][:default] = 1e-6
    end

    case_cal = case[1:n_steps]
    parameters = setup_parameters(case_cal.model)
    opt_setup = setup_parameter_optimization(
        case_cal.model, case_cal.state0, parameters, case_cal.dt, 
        case_cal.forces, objective, opt_config);

    x_opt, _, _ = Fimbul.run_optimization(opt_setup; lbfgs_args...);

    case_cal = deepcopy(opt_setup.data[:case])
    params = case_cal.parameters
    data = opt_setup.data
    params = devectorize_variables!(
        params, case_cal.model, x_opt, data[:mapper], config = data[:config])

    case_cal = JutulCase(case_cal.model, case.dt, case.forces; 
        state0 = case.state0, parameters = case_cal.parameters)

    return case_cal

end

function run_optimization(opt_setup; 
    factr = 1e8, maxfun = 200, maxiter = 200, m = 100)

    lower = opt_setup.limits.min
    upper = opt_setup.limits.max
    x0 = opt_setup.x0
    n = length(x0)
    setup = Dict(:lower => lower, :upper => upper, :x0 => copy(x0))

    prt = 1
    f! = (x) -> opt_setup.F_and_dF!(NaN, nothing, x)
    g! = (dFdx, x) -> opt_setup.F_and_dF!(NaN, dFdx, x)
    results, x_opt = LBFGSB.lbfgsb(f!, g!, x0, lb=lower, ub=upper,
        iprint = prt,
        factr = factr,
        maxfun = maxfun,
        maxiter = maxiter,
        m = m
    )

    return (x_opt, results, setup)

end