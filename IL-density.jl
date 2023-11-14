print(@__FILE__)
print("\n")
# density evaluation

# READ EXPERIMENTAL DATA
    full_file_path_rho = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-csv/rho/mix_rho.csv"]
    full_file_path_T = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-csv/rho/mix_T.csv"]

# PACKAGES
    using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
    using CSV
    using DataFrames
    using Statistics
    using PyCall
    import PyPlot

# INPUT
    p = 101.325E+03 # Pa
    exp_data = 7 # number of experimental procedures evaluated
    species1 = "C2MIMBF4"
    species2 = "C4MIMBF4"
    species3 = "C6MIMBF4"
    species4 = "C8MIMBF4"
    species5 = "C4MIMPF6"
    species6 = "C6MIMPF6"
    species7 = "C8MIMPF6"
    species8 = "C2MIMBF4_2B"
    species9 = "C4MIMBF4_2B"
    species10 = "C6MIMBF4_2B"
    species11 = "C8MIMBF4_2B"
    species12 = "C4MIMPF6_2B"
    species13 = "C6MIMPF6_2B"
    species14 = "C8MIMPF6_2B"
    species = [species1,species2,species3,species4,species5,species6,species7,species8,species9,species10,species11,species12,species13,species14]

# GENERATE MODELS
    model1 = SAFTVRMie([species1])
    model2 = SAFTVRMie([species2])
    model3 = SAFTVRMie([species3])
    model4 = SAFTVRMie([species4])
    model5 = SAFTVRMie([species5])
    model6 = SAFTVRMie([species6])
    model7 = SAFTVRMie([species7])
    model8 = SAFTVRMie([species8])
    model9 = SAFTVRMie([species9])
    model10 = SAFTVRMie([species10])
    model11 = SAFTVRMie([species11])
    model12 = SAFTVRMie([species12])
    model13 = SAFTVRMie([species13])
    model14 = SAFTVRMie([species14])
    models = [model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14]
    print(models)
    print("\n")
    n_models = length(models) # number of models evaluated
    print("n_models = ",n_models)
    print("\n")

# EXPERIMENTAL DATA
    expdata_rho = CSV.read(full_file_path_rho,DataFrame)
    expdata_T = CSV.read(full_file_path_T,DataFrame)
    print(expdata_rho)
    print("\n")
# ORGANIZE EXPERIMENTAL DATA
    exp_rho = []
    exp_T = []
    all_exp_rho = []
    all_exp_T = []
    for i ∈ 1:exp_data
        valid_rho = dropmissing(expdata_rho,[i])
        append!(exp_rho,[valid_rho[:,i]])
        append!(all_exp_rho,valid_rho[:,i])
        valid_T = dropmissing(expdata_T,[i])
        append!(exp_T,[valid_T[:,i]])
        append!(all_exp_T,valid_T[:,i])    
    end
    print("exp_rho[1] = ",exp_rho[1])
    print("\n")
# EXPERIMENTAL INTERVALS
    min_exp_rho = minimum(all_exp_rho) # g/L == kg/m³
    print("min_exp_rho = ",min_exp_rho)
    print("\n")
    max_exp_rho = maximum(all_exp_rho) # g/L == kg/m³
    print("max_exp_rho = ",max_exp_rho)
    print("\n")
    min_exp_T = minimum(all_exp_T) # K
    print("min_exp_T = ",min_exp_T)
    print("\n")
    max_exp_T = maximum(all_exp_T) # K
    print("max_exp_T = ",max_exp_T)
    print("\n")
    delta_T = max_exp_T - min_exp_T # K
    print("delta_T = ",delta_T)
    print("\n")
    length_T = trunc(Int,delta_T) # convert delta_T to integer
    print("length_T = ",length_T)
    print("\n")

#####
####
###
##
# OBTAIN critical point with 'crit_pure' function
    # crit = []
    # print("(Tc [K], Pc [Pa], Vc [m³]) = ")
    # print("\n")
    # for i ∈ 1:n_models
    #     append!(crit,[crit_pure.(models[i])])
    #     print(species[i],": ",crit[i])
    #     print("\n")
    # end

#####
####
###
##
# OBTAIN mass density with 'mass_density' function
    T = [] # K
    model_rho = [] # g/L == kg/m³
    for i ∈ 1:n_models
        append!(T,[range(min_exp_T,max_exp_T,length_T)])   
        m_d = mass_density.(models[i], p, T[i])
        append!(model_rho,[[m_d[i][1] for i ∈ 1:length_T]])
    end

# PLOT
    PyPlot.clf()
    PyPlot.rc("font", family="times new roman")
    PyPlot.figure(dpi=311)
# PLOT MODEL
    PyPlot.plot(model_rho[1],T[1],label="",linestyle="-",color="black")
    PyPlot.plot(model_rho[2],T[2],label="",linestyle="-",color="royalblue")
    PyPlot.plot(model_rho[3],T[3],label="",linestyle="-",color="hotpink")
    PyPlot.plot(model_rho[4],T[4],label="",linestyle="-",color="forestgreen")
    PyPlot.plot(model_rho[5],T[5],label="",linestyle="-",color="orange")
    PyPlot.plot(model_rho[6],T[6],label="",linestyle="-",color="indigo")
    PyPlot.plot(model_rho[7],T[7],label="",linestyle="-",color="red")
    #
    # PyPlot.plot(model_rho[1],T[1],label="",linestyle="dotted",color="black",linewidth="2")
    # PyPlot.plot(model_rho[2],T[2],label="",linestyle="dotted",color="royalblue",linewidth="2")
    # PyPlot.plot(model_rho[3],T[3],label="",linestyle="dotted",color="hotpink",linewidth="2")
    # PyPlot.plot(model_rho[4],T[4],label="",linestyle="dotted",color="forestgreen",linewidth="2")
    # PyPlot.plot(model_rho[5],T[5],label="",linestyle="dotted",color="orange",linewidth="2")
    # PyPlot.plot(model_rho[6],T[6],label="",linestyle="dotted",color="indigo",linewidth="2")
    # PyPlot.plot(model_rho[7],T[7],label="",linestyle="dotted",color="red",linewidth="2")
    #
    PyPlot.plot(model_rho[8],T[8],label="",linestyle="--",color="black")
    PyPlot.plot(model_rho[9],T[9],label="",linestyle="--",color="royalblue")
    PyPlot.plot(model_rho[10],T[10],label="",linestyle="--",color="hotpink")
    PyPlot.plot(model_rho[11],T[11],label="",linestyle="--",color="forestgreen")
    PyPlot.plot(model_rho[12],T[12],label="",linestyle="--",color="orange")
    PyPlot.plot(model_rho[13],T[13],label="",linestyle="--",color="indigo")
    PyPlot.plot(model_rho[14],T[14],label="",linestyle="--",color="red")
    # PyPlot.plot(model_rho[8],T[8],label="",linestyle="-",color="teal")
    # PyPlot.plot(model_rho[9],T[9],label="",linestyle="-",color="sienna")
    # PyPlot.plot(model_rho[10],T[10],label="",linestyle="-",color="olive")
    # PyPlot.plot(model_rho[11],T[11],label="",linestyle="-",color="dimgray")
# PLOT EXPERIMENTAL
    # https://matplotlib.org/stable/api/markers_api.html
    PyPlot.plot(exp_rho[1],exp_T[1],label="[C₂mim][BF₄]",linestyle="",marker="^",color="black")
    PyPlot.plot(exp_rho[2],exp_T[2],label="[C₄mim][BF₄]",linestyle="",marker="o",color="royalblue")
    PyPlot.plot(exp_rho[3],exp_T[3],label="[C₆mim][BF₄]",linestyle="",marker="s",color="hotpink")
    PyPlot.plot(exp_rho[4],exp_T[4],label="[C₈mim][BF₄]",linestyle="",marker="x",color="forestgreen")
    PyPlot.plot(exp_rho[5],exp_T[5],label="[C₄mim][PF₆]",linestyle="",marker="<",color="orange")
    PyPlot.plot(exp_rho[6],exp_T[6],label="[C₆mim][PF₆]",linestyle="",marker="*",color="indigo")
    PyPlot.plot(exp_rho[7],exp_T[7],label="[C₈mim][PF₆]",linestyle="",marker=">",color="red")
    # PyPlot.plot(exp_rho[8],exp_T[8],label="8",linestyle="",marker="3",color="teal")
    # PyPlot.plot(exp_rho[9],exp_T[9],label="9",linestyle="",marker="D",color="sienna")
    # PyPlot.plot(exp_rho[10],exp_T[10],label="10",linestyle="",marker="",color="olive")
    # PyPlot.plot(exp_rho[11],exp_T[11],label="11",linestyle="",marker="",color="dimgray")
# PLOT AXES
    PyPlot.legend(loc="best",frameon=false,fontsize=14)
    # PyPlot.legend(loc="best",frameon=false,fontsize=12,bbox_to_anchor=(0.5,0.5))
    PyPlot.xlabel("Density (kg/m³)",fontsize=16)
    PyPlot.ylabel("Temperature (K)",fontsize=16)
    PyPlot.xticks(fontsize=12)
    PyPlot.yticks(fontsize=12)
    # PyPlot.xlim([min_exp_rho,max_exp_rho])
    PyPlot.xlim([min_exp_rho,1520])
    PyPlot.ylim([min_exp_T,max_exp_T])
display(PyPlot.gcf())
