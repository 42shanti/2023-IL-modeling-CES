print(@__FILE__)
print("\n")
# BUB point AAD
# calculate bub point in the same composition of experimental data

# PACKAGES
    using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
    using CSV
    using DataFrames
    using Statistics
    using PyCall
    import PyPlot

# only 'INPUT' and 'EXPERIMENTAL PATH' need to be modified

# INPUT
    species1 = "CO2"
    # species2 = "C2MIMSCN"
    species2 = "C2MIMSCN_2B"
    T_ = 363.15 # K
    # exp_column = 1 # 'exp_column = 1' refers to experimental data of T1
    exp_column = 2 # 'exp_column = 2' refers to experimental data of T2
    # x1 = [0.169, 0.251, 0.318, 0.364, 0.410] # test
    # X1 = Clapeyron.FractionVector.(x1) # test
    # print("X = ",X1," #length = ",length(X1)) # test
    # print("\n")

# GENERATE MODELS
    model1 = SAFTVRMie([species1,species2])
    print(model1)
    print("\n")

# EXPERIMENTAL PATH
    full_file_path_x = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-csv/bub-co2/x-C2MIMSCN.csv"] # composition file
    full_file_path_p = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-csv/bub-co2/p-C2MIMSCN.csv"] # pressure file
# READ experimental data
    expdata_x = CSV.read(full_file_path_x,DataFrame)
    expdata_p = CSV.read(full_file_path_p,DataFrame)
    print(expdata_x)
    print("\n\n")
    print("collumn 1 = pressure in Temperature 1; collumn 2 = pressure in T2:")
    print("\n")
    print(expdata_p)
    print("\n")
# ORGANIZE experimental data
    exp_x = [] # mole fraction
    exp_p = [] # Pa
    all_exp_x = []
    all_exp_p = []
    n_columns = 2 # number of collumns. T1 = collumn 1; T2 = column 2
    # the following 'for' is to remove missing values to make consistent plot
    for i ∈ 1:n_columns
        valid_x = dropmissing(expdata_x,[i]) # remove missing value
        append!(exp_x,[valid_x[:,i]])
        append!(all_exp_x,valid_x[:,i])
        valid_p = dropmissing(expdata_p,[i])
        append!(exp_p,[valid_p[:,i]])
        append!(all_exp_p,valid_p[:,i])
    end
# BOUNDARIES
    # pressure
    min_exp_p = minimum(all_exp_p) # Pa
    # print("min_exp_p = ",min_exp_p)
    # print("\n")
    max_exp_p = maximum(all_exp_p) # Pa
    # print("max_exp_p = ",max_exp_p)
    # print("\n")
    # mole fraction
    min_exp_x = minimum(all_exp_x) # mole fraction
    print("min_exp_x = ",min_exp_x)
    print("\n")
    max_exp_x = maximum(all_exp_x) # mole fraction
    print("max_exp_x = ",max_exp_x)
    print("\n")
# ORGANIZE BOUNDARIES
    p_begin = min_exp_p*1E-05
    print("p_begin = ",p_begin)
    print("\n")
    p_end = max_exp_p*1E-05
    print("p_end = ",p_end)
    print("\n")

    x1 = exp_x[exp_column]
    X1 = Clapeyron.FractionVector.(x1)
    print("X = ",X1," #length = ",length(X1))
    print("\n")

# SIMPLIFIED BUBBLE PRESURE MODEL
    bub1 = []
    p1 = []
    for i ∈ 1:length(x1)
        bp = bubble_pressure(model1,T_,X1[i])
        append!(p1,bp[1])
    end
    print(p1)
    print("\n")

# AAD
model_median = median(p1)
print("model_median = ",model_median)
print("\n")
exp_median = median(exp_p[exp_column])
print("exp_median = ",exp_median)
print("\n")
AAD = 100*((model_median-exp_median)/exp_median)
print("%AAD = ",abs(AAD))
print("\n")

# PLOT
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(dpi=311)
# plot model
    # T_
    PyPlot.plot(x1,p1*1E-05,label="",linestyle="solid",marker="*",color="black")
# plot experimental
    # T_
    PyPlot.plot(exp_x[exp_column],exp_p[exp_column]*1E-05,label=string(T_," ","K"),linestyle="",marker="o",color="black")
# plot axes
    PyPlot.legend(loc="best",frameon=false,fontsize=16) 
    PyPlot.xlabel("Mole fraction",fontsize=16)
    PyPlot.ylabel("Pressure (bar)",fontsize=16)
    PyPlot.xticks(fontsize=12)
    PyPlot.yticks(fontsize=12)
# limits
    # PyPlot.xlim([min_exp_x,max_exp_x])
    # PyPlot.xlim([x_begin,x_end])
    # PyPlot.ylim([p_begin,p_end])
    # PyPlot.xlim([y1[1][begin],y1[1][end]])
    # PyPlot.ylim([p1[1][begin]*1E-05,p1[1][end]]*1E-05)
    # PyPlot.xlim([0.0,1.1*max_exp_x])
    # PyPlot.ylim([3.0,1.1*p_end])
    # PyPlot.xlim([0.0,0.4])
    # PyPlot.ylim([3.0,151])
#
display(PyPlot.gcf())
