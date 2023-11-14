print(@__FILE__)
print("\n")
# @time begin
# Parameter Estimation (PE)
# for the full PE explanation follow the source:
# https://github.com/ClapeyronThermo/Clapeyron.jl/tree/master/examples

# EXPERIMENTAL PATHs CSV
    full_file_path_rho = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-redux/csv-pe-rho/C4MIMBF4.csv"

# PACKAGES
    using Clapeyron # https://doi.org/10.1021/acs.iecr.2c00326
    using Statistics
    using BlackBoxOptim

# GENERATE MODEL you wish to parametrize
    species = "C4MIMBF4"
    model = SAFTVRMie([species]);
    print(model)
    print("\n")

# ESTIMATION parameters
    segmentlower = 5.42
    segmentupper = 5.52
    segmentguess = median([segmentlower,segmentupper])
    sigmalower = 2.0
    sigmaupper = 5.0
    sigmaguess = median([sigmalower,sigmaupper])
    epsilonlower = 200.0
    epsilonupper = 500.0
    epsilonguess = median([epsilonlower,epsilonupper])
    # epsilon_assoclower = 1000.0
    # epsilon_assocupper = 3000.0
    # epsilon_assocguess = median([epsilon_assoclower,epsilon_assocupper])
# ESTIMATION framework
    toestimate = [
        Dict(
            :param => :segment,
            # :indices => 1,
            # :cross_assoc => true,
            :lower => segmentlower,
            :upper => segmentupper,
            :guess => segmentguess
        ),
        Dict(
            :param => :sigma, # [Å]
            # :indices => 1,
            # :cross_assoc => true,
            :factor => 1E-10, # convert [Å] to SI unit [m]
            :lower => sigmalower,
            :upper => sigmaupper,
            :guess => sigmaguess
        ),
        Dict(
            :param => :epsilon, # [K]
            # :indices => 1,
            # :cross_assoc => true,
            :lower => epsilonlower,
            :upper => epsilonupper,
            :guess => epsilonguess
        ),
        # Dict(
        #     :param => :epsilon_assoc,
        #     :lower => epsilon_assoclower,
        #     :upper => epsilon_assocupper,
        #     :guess => epsilon_assocguess
        # ),
    ];

# PROPERTY function you with to estimate
    # for more information follow https://doi.org/10.1021/acs.iecr.2c00326
    function mass_rho(model::EoSModel,p,T)
        md = mass_density(model,p,T)
        return md[1]
    end

# ESTIMATOR function
    estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_rho]);

# OPTMIZATION
    nparams = length(initial)
    bounds  = [(lower[i],upper[i]) for i in 1:nparams]

    result = BlackBoxOptim.bboptimize(objective; 
            SearchRange = bounds, 
            NumDimensions = nparams,
            MaxSteps=10000,
            PopulationSize = 1000,
            TraceMode=:silent)

    params = BlackBoxOptim.best_candidate(result);

# PRINT PARAMETERS
    print("[segment, sigma, epsilon]=",params)
    print("\n")

# end#@time
