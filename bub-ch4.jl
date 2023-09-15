print(@__FILE__)
print("\n")
@time begin
# bub aad evaluation
# main code - here we generate the solubility graphs

# PACKAGES
    using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
    using CSV
    using DataFrames
    using Statistics
    using PyCall
    import PyPlot

# EXPERIMENTAL data path
    # full_file_path_x = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03/csv-aad-bub/x-C8MIMTF2N.csv"]
    # full_file_path_p = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03/csv-aad-bub/p-C8MIMTF2N.csv"]
    # full_file_path_x = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03 forestgreenux/csv-aad-bub/x-C8MIMTF2N.csv"]
    # full_file_path_p = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03 forestgreenux/csv-aad-bub/p-C8MIMTF2N.csv"]
    full_file_path_x = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03 csv/bub-ch4/x-C8MIMTF2N.csv"]
    full_file_path_p = ["C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp brainstorm/brainstorm 03 csv/bub-ch4/p-C8MIMTF2N.csv"]

# INPUT
    x_length = 33 # interval points - impacts runtime and convergence
    T1 = 313.13 # K
    T2 = 353.15 # K
    species1 = "CH4"
    species2 = "C8MIMTF2N"
    species3 = "C8MIMTF2N2B"

# GENERATE models
    model1 = SAFTVRMie([species1,species2])
    model2 = SAFTVRMie([species1,species3])
    models = [model1,model2];
    # models = [1]
    print(models)
    print("\n")
    n_models = length(models) # number of models evaluated
    print("n_models = ",n_models)
    print("\n")

# READ experimental data
    expdata_x = CSV.read(full_file_path_x,DataFrame)
    expdata_p = CSV.read(full_file_path_p,DataFrame)
    print(expdata_x)
    print("\n")
    print(expdata_p)
    print("\n")
    # T = [303.15] # K
    # T1 = T[1]
    T = [T1 T2] # K
    length_T = length(T)
    # print("length_T = ",length_T)
    # print("\n")
# ORGANIZE experimental data
    exp_x = [] # mole fraction
    exp_p = [] # Pa
    all_exp_x = []
    all_exp_p = []
    for i ∈ 1:length_T
        valid_x = dropmissing(expdata_x,[i]) # remove missing value
        append!(exp_x,[valid_x[:,i]])
        append!(all_exp_x,valid_x[:,i])
        valid_p = dropmissing(expdata_p,[i])
        append!(exp_p,[valid_p[:,i]])
        append!(all_exp_p,valid_p[:,i])
        # print(T[i]," K - exp_x = ",exp_x[i])
        # print("\n")
        # print(T[i]," K - exp_p = ",exp_p[i])
        # print("\n") 
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
# ORGANIZE boundaries
    p_begin = min_exp_p*1E-05
    print("p_begin = ",p_begin)
    print("\n")
    p_end = max_exp_p*1E-05
    print("p_end = ",p_end)
    print("\n")
    #
    x_begin = 0.0
    x_end = 1.0
    # x_begin = min_exp_x
    # x_end = max_exp_x
    ###
    ##
    #
    # x_length = 17
    x = range(x_begin,x_end,x_length)
    X = Clapeyron.FractionVector.(x)
    print("X = ",X," #length = ",length(X))
    print("\n")

#####
####
###
##
# MODEL the 'pxy' envelope using the 'bubble_pressure' function
    # T1 = 309.15 # K
    # --- T1 ---
    y1 = []
    p1 = []
    for i ∈ 1:n_models
        bub1 = []
        v0 =[]
        for j ∈ 1:x_length
            if j==1
                append!(bub1, [bubble_pressure(models[i],T1,X[j])])
                v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
            else
                append!(bub1, [bubble_pressure(models[i],T1,X[j];v0=v0)])
                v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
            end
        end
        append!(y1,[append!([bub1[i][4][1] for i ∈ 1:x_length],reverse(x))])
        append!(p1,[append!([bub1[i][1] for i ∈ 1:x_length],[bub1[i][1] for i ∈ x_length:-1:1])])
    end
    ###
    ##
    # --- T2 ---
        # x_length = 30
        # x = range(x_begin,x_end,x_length)
        # X = Clapeyron.FractionVector.(x)
        # print("X = ",X," #length = ",length(X))
        # print("\n")
    ###
    ##
    #
    y2 = []
    p2 = []
    for i ∈ 1:n_models
        bub2 = []
        v0 = []
        for j ∈ 1:x_length
            if j==1
                append!(bub2, [bubble_pressure(models[i],T2,X[j])])
                v0 = [log10(bub2[j][2]),log10(bub2[j][3]),bub2[j][4][1],bub2[j][4][2]]
            else
                append!(bub2, [bubble_pressure(models[i],T2,X[j];v0=v0)])
                v0 = [log10(bub2[j][2]),log10(bub2[j][3]),bub2[j][4][1],bub2[j][4][2]]
            end
        end
        append!(y2,[append!([bub2[i][4][1] for i ∈ 1:x_length],reverse(x))])
        append!(p2,[append!([bub2[i][1] for i ∈ 1:x_length],[bub2[i][1] for i ∈ x_length:-1:1])])
    end

###
##
# PRINT RESULTS - saving the terminal results is advised  
    print("--- Model results ---")
    print("\n")
    print("species: ",species1,"-",species2)
    print("\n")
    print("T1 = ",T1)
    print("\n")
    print("T2 = ",T2)
    print("\n")
    print("p1[1] = ",p1[1]," #length = ",length(p1[1]))
    print("\n")
    print("p1[2] = ",p1[2]," #length = ",length(p1[2]))
    print("\n")
    print("y1[1] = ",y1[1])
    print("\n")
    print("y1[2] = ",y1[2])
    print("\n")
    print("p2[1] = ",p2[1]," #length = ",length(p2[1]))
    print("\n")
    print("p2[2] = ",p2[2]," #length = ",length(p2[2]))
    print("\n")
    print("y2[1] = ",y2[1])
    print("\n")
    print("y2[2] = ",y2[2])
    print("\n")
    # print("p3[1] = ",p3[1])
    # print("\n")
    # print("y3[1] = ",y3[1])
    # print("\n")

###
##
# PLOT
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
# plot model
    # T1
    PyPlot.plot(y1[1],p1[1]*1E-05,label="",linestyle="solid",color="black")
    PyPlot.plot(y1[2],p1[2]*1E-05,label="",linestyle="dashed",color="black")
    # PyPlot.plot(y3[1][35:53],p3[1][35:53]*1E-05,label="",linestyle="solid",color="black")
    # T2
    PyPlot.plot(y2[1],p2[1]*1E-05,label="",linestyle="-",color="forestgreen")
    PyPlot.plot(y2[2],p2[2]*1E-05,label="",linestyle="--",color="forestgreen")
# plot experimental
    # T1
    PyPlot.plot(exp_x[1],exp_p[1]*1E-05,label="313.15 K",linestyle="",marker="o",color="black")
    # T2
    PyPlot.plot(exp_x[2],exp_p[2]*1E-05,label="353.15 K",linestyle="",marker="v",color="forestgreen")
# plot axes
    PyPlot.legend(loc="best",frameon=false,fontsize=16) 
    PyPlot.xlabel("CH₄ mole fraction",fontsize=16)
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
    # PyPlot.ylim([0.0,1.1*p_end])
    PyPlot.xlim([0.0,0.179])
    PyPlot.ylim([4.0,63])
#
display(PyPlot.gcf())

# species2 = "C8MIMTF2N"

end#@time

#     x_length = 27
#     x = range(x_begin,x_end,x_length)
#     X = Clapeyron.FractionVector.(x)
#     print("X = ",X," #length = ",length(X))
#     print("\n")
# # MODEL the 'pxy' envelope using the 'bubble_pressure' function
#     T1 = 299.85 # K
#     # --- T1 ---
#     y3 = []
#     p3 = []
#     for i ∈ 1:n_models
#         bub1 = []
#         v0 =[]
#         for j ∈ 1:x_length
#             if j==1
#                 append!(bub1, [bubble_pressure(model1,T1,X[j])])
#                 v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
#             else
#                 append!(bub1, [bubble_pressure(model1,T1,X[j];v0=v0)])
#                 v0 = [log10(bub1[j][2]),log10(bub1[j][3]),bub1[j][4][1],bub1[j][4][2]]
#             end
#         end
#         append!(y3,[append!([bub1[i][4][1] for i ∈ 1:x_length],reverse(x))])
#         append!(p3,[append!([bub1[i][1] for i ∈ 1:x_length],[bub1[i][1] for i ∈ x_length:-1:1])])
#     end

