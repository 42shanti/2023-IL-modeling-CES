print(@__FILE__)
print("\n")
# code to find the constants of the molecular volume and SAFT correlation

# PACKAGES
    using CSV
    using DataFrames
    using Plots
    using LaTeXStrings
    using PyCall
    import PyPlot

# DATA PATH
    full_file_path_2B = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-redux/correlation-2B.csv"
    full_file_path_4C = "C:/Users/bc_usp/My Drive (cleitonberaldo@usp.br)/usp-brainstorm/brainstorm-03-redux/correlation-4C.csv"

# READ DATA
    correlation_2B = CSV.read([full_file_path_2B],DataFrame)
    correlation_4C = CSV.read([full_file_path_4C],DataFrame)
    print(correlation_2B)
    print("\n")
    print(correlation_4C)
    print("\n")

    c1_lim_2B = [correlation_2B[begin,1],correlation_2B[end,1]]
    c2_lim_2B = [correlation_2B[begin,2],correlation_2B[end,2]]
    c3_lim_2B = [correlation_2B[begin,3],correlation_2B[end,3]]
    c4_lim_2B = [correlation_2B[begin,4],correlation_2B[end,4]]

    c1_lim_4C = [correlation_4C[begin,1],correlation_4C[end,1]]
    c2_lim_4C = [correlation_4C[begin,2],correlation_4C[end,2]]
    c3_lim_4C = [correlation_4C[begin,3],correlation_4C[end,3]]
    c4_lim_4C = [correlation_4C[begin,4],correlation_4C[end,4]]

# PLOT
    PyPlot.clf()
    plot_font = "times new roman"
    PyPlot.rc("font", family=plot_font)
    PyPlot.figure(dpi=311)

# PLOT 1 - m (dimensionless)

    # 2B
    # PyPlot.plot(correlation_2B[:,1],correlation_2B[:,2],label="m = 0.00987×V"*L"\rm{_m}"*" + 3.1266",linestyle="",marker="o",color="black")
    # PyPlot.plot(c1_lim_2B,c2_lim_2B,label="",linestyle="--",marker="",color="black")
    
    # 4C
    PyPlot.plot(correlation_4C[:,1],correlation_4C[:,2],label="m = 0.00944×V"*L"\rm{_m}"*" + 4.8950",linestyle="",marker="o",color="black")
    PyPlot.plot(c1_lim_4C,c2_lim_4C,label="",linestyle="-",marker="",color="black")

# PLOT 2 - m*σ³ (Å³)

    # 2B
    # PyPlot.plot(correlation_2B[:,1],correlation_2B[:,3]./100,label="m σ³ = 1.1336×V"*L"\rm{_m}"*" – 23.281",linestyle="",marker="v",color="forestgreen")
    # PyPlot.plot(c1_lim_2B,c3_lim_2B./100,label="",linestyle="--",marker="",color="forestgreen")

    # 4C
    PyPlot.plot(correlation_4C[:,1],correlation_4C[:,3]./100,label="m σ³ = 1.1396×V"*L"\rm{_m}"*" – 28.170",linestyle="",marker="v",color="forestgreen")
    PyPlot.plot(c1_lim_4C,c3_lim_4C./100,label="",linestyle="-",marker="",color="forestgreen")

# PLOT 3 - m*ϵ (K)

    # 2B
    # PyPlot.plot(correlation_2B[:,1],correlation_2B[:,4]./1000,label="m ϵ = 3.4073×V" * L"\rm{_m}" * " + 1101.95",linestyle="",marker="^",color="hotpink")
    # PyPlot.plot(c1_lim_2B,c4_lim_2B./1000,label="",linestyle="--",marker="",color="hotpink")

    # 4C
    PyPlot.plot(correlation_4C[:,1],correlation_4C[:,4]./1000,label="m ϵ = 3.8892×V" * L"\rm{_m}" * " + 1195.56",linestyle="",marker="^",color="hotpink")
    PyPlot.plot(c1_lim_4C,c4_lim_4C./1000,label="",linestyle="-",marker="",color="hotpink")

# PLOT AXES
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html
    PyPlot.legend(loc=6,frameon=false,fontsize=14)
    PyPlot.xlabel("Molecular volume (Å³)",fontsize=16)
    PyPlot.ylabel("m, m×σ³ (Å³)/100, m×ϵ (K)/1000",fontsize=16)
    PyPlot.xticks(fontsize=12)
    PyPlot.yticks(fontsize=12)

display(PyPlot.gcf())
