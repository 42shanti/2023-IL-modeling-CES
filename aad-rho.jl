print(@__FILE__)
print("\n")
# density AAD calculation

# PACKAGES
    using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
    using CSV
    using DataFrames
    using Statistics
    using PyCall
    import PyPlot

# INPUT
    species1 = "C2MIMBF4"
    p = 101.325E+03 # Pa

# GENERATE MODELS
    model1 = SAFTVRMie([species1])
    models = [model1]
    print(models)
    print("\n")
    n_models = length(models) # number of models evaluated
    print("n_models = ",n_models)
    print("\n")

    
# EXPERIMENTAL DATA PATH
    # full_file_path = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03/csv-aad-rho/C2MIMBF4.csv"]
    full_file_path = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-redux/csv-aad-rho/C2MIMBF4.csv"]

# EXPERIMENTAL DATA
    expdata = CSV.read(full_file_path,DataFrame)
    print(expdata)
    print("\n")
    exp_rho1st = expdata[:,1] # all values from collumn 1
# EXPERIMENTAL INTERVALS
    T_begin = expdata[begin,2] # K
    print("T_begin = ",T_begin)
    print("\n")
    T_end = expdata[end,2] # K
    print("T_end = ",T_end)
    print("\n")
    delta_T = T_end - T_begin # K
    print("delta_T = ",delta_T)
    print("\n")
    length_T = trunc(Int,delta_T) # convert delta_T to integer
    print("length_T = ",length_T)
    print("\n")
    rho_begin = expdata[begin,1] # g/L == kg/m³
    print("rho_begin = ",rho_begin)
    print("\n")
    rho_end = expdata[end,1] # g/L == kg/m³
    print("rho_end = ",rho_end)
    print("\n")

# OBTAIN mass density with 'mass_density' method
    T = [] # K
    model_rho = [] # g/L == kg/m³
    for i ∈ 1:n_models
        append!(T,[range(T_begin,T_end,length_T)])   
        m_d = mass_density.(models[i], p, T[i])
        append!(model_rho,[[m_d[i][1] for i ∈ 1:length_T]])
    end

# AAD
    model_rho_median = median(model_rho[1])
    print("model_rho_median = ",model_rho_median)
    print("\n")
    exp_rho_median = median(expdata[!,1])
    print("exp_rho_median = ",exp_rho_median)
    print("\n")
    AAD = 100*((model_rho_median-exp_rho_median)/exp_rho_median)
    print("#######")
    print("\n")
    print("AAD% = ",abs(AAD))
    print("\n")

# PLOT
    PyPlot.clf()
    PyPlot.rc("font", family="times new roman")
# PLOT MODEL
    PyPlot.plot(model_rho[1],T[1],label="SAFT",linestyle="-",color="black")
# PLOT EXPERIMENTAL
    PyPlot.plot(expdata[!,1],expdata[!,2],label="Experimental",linestyle="",marker="o",color="black")
# PLOT AXES
    PyPlot.legend(loc="best",frameon=false,fontsize=12)
    # PyPlot.legend(loc="best",frameon=false,fontsize=12,bbox_to_anchor=(0.5,0.5))
    PyPlot.xlabel("Density (kg/m³)",fontsize=16)
    PyPlot.ylabel("Temperature (K)",fontsize=16)
    PyPlot.xticks(fontsize=12)
    PyPlot.yticks(fontsize=12)
    # PyPlot.xlim([rho_end,rho_begin])
    # PyPlot.ylim([T_begin,T_end])
    PyPlot.xlim([model_rho[1][end],model_rho[1][begin]])
    PyPlot.ylim([T[1][begin],T[1][end]])
display(PyPlot.gcf())
