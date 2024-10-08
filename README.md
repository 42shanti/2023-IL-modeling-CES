# 2023 Scripts - IL solubility modeling - Chemical Engineering Science
Paper: https://doi.org/10.1016/j.ces.2023.119610

# Codes used to model the IL systems:


IL-parameter-estimation.jl - parameter estimation

IL-density.jl - obtain mass density

IL-bub-pressure.jl - VLE calculation (old version - does not work on newer versions of Clapeyron.jl)

IL-bub-pressure-2024-05.jl - VLE calculation (updated script - works just fine in Clapeyron.jl version 0.6.2)


# Julia packages version
name = "Clapeyron"
uuid = "7c7805af-46cc-48c9-995b-ed0ed2dc909a"
authors = ["Hon Wa Yew <yewhonwa@gmail.com>", "Pierre Walker <pwalker@mit.edu>", "Andrés Riedemann <andres.riedemann@gmail.com>"]
version = "0.4.0"

[deps]
BlackBoxOptim = "a134a8b2-14d6-55f6-9291-3336d3ab0209"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
PackedVectorsOfVectors = "7713531c-48ef-4bdd-9821-9ff7a8736089"
PositiveFactorizations = "85a6dd25-e78a-55b7-8502-1745935b8125"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Scratch = "6c6a2e73-6563-6170-7368-637461726353"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
ThermoState = "e7b6519d-fdf7-4a33-b544-2b37a9c1234a"
UUIDs = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
BlackBoxOptim = "^0.6.2"
CSV = "0.10"
DiffResults = "^1.1"
FillArrays = "^0.12, 0.13"
ForwardDiff = "^0.10"
LogExpFunctions = "^0.2,^0.3"
NLSolvers = "0.4"
PackedVectorsOfVectors = "^0.1.2"
PositiveFactorizations = "^0.2"
Roots = "^2.0.4"
Scratch = "^1.1"
StaticArrays = "^1.5.9"
Tables = "^1.8"
ThermoState = "^0.5"
Unitful = "^1.12"
julia = "1.6"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
