print(@__FILE__)
print("\n")
@time begin
# Parameter Estimation (PE)
# for the full PE explanation follow the source:
# https://github.com/ClapeyronThermo/Clapeyron.jl/tree/master/examples

# EXPERIMENTAL PATHs CSV
    full_file_path_rho = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-redux/csv-pe-rho/C8MIMPF6.csv"
    # full_file_path_sound = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03/csv-pe-sound/C8MIMPF6.csv"
    # full_file_path_cp = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03/csv-pe-cp/C8MIMPF6.csv"
    # full_file_path_mrho = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp clapeyron/Clapeyron.jl/2023-05-brainstorm-03/1st-test/pe-multi-rho.csv"
    # full_file_path_sat = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03/csv-pe-sat/C8MIMPF6.csv"

# PACKAGES
    using Clapeyron # https://doi.org/10.1021/acs.iecr.2c00326
    using Statistics
    using BlackBoxOptim

# GENERATE MODEL you wish to parametrize
    model = SAFTVRMie(["C8MIMPF6"]);
    # model = SAFTVRMie(["co2","C8MIMPF6"]);
    print(model)
    print("\n")

# ESTIMATION parameters
    segmentlower = 6.92
    segmentupper = 7.02
    segmentguess = median([segmentlower,segmentupper])
    sigmalower = 2.0
    sigmaupper = 5.0
    sigmaguess = median([sigmalower,sigmaupper])
    epsilonlower = 200.0
    epsilonupper = 500.0
    epsilonguess = median([epsilonlower,epsilonupper])
    epsilon_assoclower = 1000.0
    epsilon_assocupper = 3000.0
    epsilon_assocguess = median([epsilon_assoclower,epsilon_assocupper])
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
    # function sat_p(model::EoSModel,T)
    #     sat = saturation_pressure(model,T)
    #     return sat[1]
    # end
    # function sound(model::EoSModel,p,T)
    #     ss = speed_of_sound(model,p,T)
    #     return ss[2]
    # end
    # function heat_capacity(model::EoSModel,p,T)
    #     cp = isobaric_heat_capacity(model,p,T)
    #     return cp[1]
    # end
    # function mass_rho(model::EoSModel,p,T)
    #     md = mass_density(model,p,T)
    #     return 1/md[2]
    # end

print("function ok")
print("\n")

# ESTIMATOR function
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[(2.,"C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp clapeyron/Clapeyron.jl/worked-examples/data/saturation_pressure.csv"),(1.,"C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp clapeyron/Clapeyron.jl/worked-examples/data/saturation_liquid_density.csv")]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[(99.,full_file_path_rho),(1.,full_file_path_sound)]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_sound]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_cp]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_mrho]);
    estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_rho]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[(999.,full_file_path_rho),(1.,full_file_path_sat)]);
    # estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_sat]);

print("estimator ok")
print("\n")

nparams = length(initial)
bounds  = [(lower[i],upper[i]) for i in 1:nparams]

result = BlackBoxOptim.bboptimize(objective; 
        SearchRange = bounds, 
        NumDimensions = nparams,
        MaxSteps=10000,
        PopulationSize = 1000,
        TraceMode=:silent)

params = BlackBoxOptim.best_candidate(result);

print("[segment, sigma, epsilon]=",params)
# print("[segment, sigma, epsilon, epsilon_assoc]=",params)
# print("[epsilon, epsilon_assoc]=",params)
print("\n")

end#@time
