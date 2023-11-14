print(@__FILE__)
print("\n")
@time begin
# bub evaluation
# main code - here we calcualte the solubility

# PACKAGES
    using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
    using CSV
    using DataFrames
    using Statistics
    using PyCall
    import PyPlot

# EXPERIMENTAL DATA PATH
    full_file_path_x = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/bs03-r3/bs03r3-csv/x-C2MIMSCN.csv"] # experimental fraction
    full_file_path_p = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/bs03-r3/bs03r3-csv/p-C2MIMSCN.csv"] # experimental pressure

# INPUT
    x_length_1 = 12 # interval points - impacts runtime and convergence
    # numerical pitftalls can be avoided with a well defined x_length
    T1 = 313.15
    experimental_datasets = 3 # number of experimental datasets
    # ie, three temperatures have 3 fractions and 3 pressures
    species0 = "CO2"
    species1 = "C2MIMSCN"

# GENERATE MODELS
    model1 = SAFTVRMie([species0,species1])
    models = [model1];
    print(models)
    print("\n")
    n_models = length(models) # number of models evaluated

# READ EXPERIMENTAL DATA
    expdata_x = CSV.read(full_file_path_x,DataFrame)
    expdata_p = CSV.read(full_file_path_p,DataFrame)
    # print(expdata_x)
    # print("\n")
    # print(expdata_p)
    # print("\n")
    T = [T1] # K
# ORGANIZE EXPERIMENTAL DATA
    exp_x = [] # mole fraction
    exp_p = [] # kPa
    all_exp_x = []
    all_exp_p = []
    # column 1 = experimental data set of temperature 1; column 2 = exp data set of temperature 2...
    for i ∈ 1:experimental_datasets
        valid_x = dropmissing(expdata_x,[i]) # remove missing values
        append!(exp_x,[valid_x[:,i]])
        append!(all_exp_x,valid_x[:,i])
        valid_p = dropmissing(expdata_p,[i])
        append!(exp_p,[valid_p[:,i]])
        append!(all_exp_p,valid_p[:,i])
    end

# BOUNDARIES
    # PRESSURE
    min_exp_p = minimum(all_exp_p) # kPa
    max_exp_p = maximum(all_exp_p) # kPa
    min_exp_x = minimum(all_exp_x) # mole fraction
    max_exp_x = maximum(all_exp_x) # mole fraction
# ORGANIZE BOUNDARIES
    p_begin = min_exp_p*1E-02 # convert pressure from kPa to bar
    p_end = max_exp_p*1E-02 # bar
    # COMPOSITION INTERVAL CALCULATION
        # in doubt leave it 0.0–1.0
    x_begin = 0.0 # mole fraction
    x_end = 1.0 # mole fraction

#
##
####
#####
######
#######################
######
#####
####
###
##
# MODEL THE 'PXY' ENVELOPE USING THE 'bubble_pressure' FUNCTION

    print("BUB PRESSURE EVALUATION")
    print("\n")
    x_1 = range(x_begin,x_end,x_length_1)
    X_1 = Clapeyron.FractionVector.(x_1) # it creates a vector for VLE calculation
    y1 = [] # store vapor fraction here
    p1 = [] # store pressure fraction here
    for i ∈ 1:n_models # evaluate all models
        bub1 = [] # store results of 'bubble_pressure' function here
        # function information: [Table 7] of https://doi.org/10.1021/acs.iecr.2c00326
        v0 =[] # magic
        for j ∈ 1:x_length_1
            if j==1
                append!(bub1, [bubble_pressure(model1,T1,X_1[j])])
                v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
            else
                append!(bub1, [bubble_pressure(model1,T1,X_1[j];v0=v0)])
                v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
            end
        end
        append!(y1,[append!([bub1[i][4][1] for i ∈ 1:x_length_1],reverse(x_1))])
        append!(p1,[append!([bub1[i][1] for i ∈ 1:x_length_1],[bub1[i][1] for i ∈ x_length_1:-1:1])])
    end
    ###
    ####
    #####
    ######
    #####
    ####
    ###
#######
######
#####
####
###
##
# PRINT RESULTS - saving the terminal results is highly advised  
    print("--- Model results ---")
    print("\n")
    print("species: ",species0,"-",species1)
    print("\n")
    print("T1 = ",T1)
    print("\n")
    print("p1 (4C-T1) = ",p1[1])
    print("\n")
    print("y1 = ",y1[1])
    print("\n")

end#@time

###########
##########
#########
########
#######
######
#####
####
###
##
# PLOT
    PyPlot.clf()
    PyPlot.rc("font", family="times new roman")
    PyPlot.figure(dpi=311)

# PLOT MODEL
    # 1 (4C-T1)
    PyPlot.plot(y1[1],p1[1].*1E-05,label="",linestyle="solid",color="black")
# PLOT EXPERIMENTAL
    # T1
    PyPlot.plot(exp_x[1],exp_p[1]*1E-02,label="313.15 K",linestyle="",marker="^",color="black")
    # T2
    # PyPlot.plot(exp_x[2],exp_p[2]*1E-02,label="343.15 K",linestyle="",marker="o",color="royalblue")
    # T3
    # PyPlot.plot(exp_x[3],exp_p[3]*1E-02,label="373.15 K",linestyle="",marker="s",color="hotpink")
# PLOT AXES
    PyPlot.legend(loc="upper left",frameon=false,fontsize=16) 
    PyPlot.xlabel("CO₂ mole fraction",fontsize=16)
    PyPlot.ylabel("Pressure (bar)",fontsize=16)
    PyPlot.xticks(fontsize=12)
    PyPlot.yticks(fontsize=12)
        # LIMITS
        x_final = max_exp_x
        P_final = p_end
        # x_final = 0.44
        # P_final = 550
        PyPlot.xlim([0.0,x_final])
        PyPlot.ylim([9.9,P_final])
        # PyPlot.text(0.07.*x_final,P_final./(1.8), "CO₂–[C₂mim][SCN]", fontsize=16)
    #
    display(PyPlot.gcf())

# species1 = "C2MIMSCN" # quick ctrl h
#
##
###₁₂₃₄₅₆₈
##
#