## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh & Affan

module RSV

## Loading necessary packages
using DataFrames
using Distributions
using StatsBase
using StaticArrays
using Random

# construct agent's types as an object
mutable struct Human
    idx::Int64
    health::Int64       # 0 = susc, 1 = infected
    age::Int64          # in yeras
    agegroup::Int64    # G1 = < 1, G2= 2, G3= 3, G4= 4, G5= 5:19, G6= > 19
    smoker::Int64       # 1 = no,  2 = yes
    preterm:: Int64     # 1 = no,  2 = yes
    Human() = new()
end

mutable struct Dwelling
    idx::Int64
    humans::Array{Int64}
    size::Int64
    Dwelling() = new()
end

######## global variables
# Nunavik age population. Downloaded from Census Canada 2016
# data organized in file 'parameters_Shokoofeh'
const agedist =  Categorical(@SVector[0.022727273, 0.026136364, 0.026136364, 0.026136364, 0.022727273, 0.310606061, 0.565530303]) # categorical probability distribution of discrete age_groups
const agebraks = @SVector[0:0.99, 1:1.99, 2:2.99, 3:3.99, 4:4.99, 5:18.99, 19:100] #age_groups
const predist = Categorical(@SVector[0.972027972, 0.027972028])  # distribution of preterm infants
const smokdist = Categorical(@SVector[0.37, 0.63])  # distribution of preterm infants
const dwellsize_dist = Categorical(@SVector[0.208275862,0.164137931,0.153103448,0.16137931,0.313103448]) # household's size dist.: 1, 2, 3, 4, 5 >
const dwellbraks = @SVector[1:1, 2:2, 3:3, 4:4, 5:15] # house size
#const P = ModelParameters() # I'll make the parameter file
const humans = Array{Human}(undef, 13200) # 13200 is total population of Nunavik
const dwells = Array{Dwelling}(undef, 1000)


export humans, dwells


############################ setting up agents for model
function init_humans()
    @inbounds for i = 1:length(humans)
        humans[i] = Human()              ## create an empty human
        humans[i].idx = i
        humans[i].health = 0
        humans[i].age = 0
        humans[i].agegroup = 0
        humans[i].smoker = 1
        humans[i].preterm = 1
    end
end
export init_humans

function init_dwellings()
    @inbounds for i = 1:length(dwells)
        dwells[i] = Dwelling()             ## create an empty house
        dwells[i].idx = i
        dwells[i].size = rand(dwellbraks[rand(dwellsize_dist)])
        dwells[i].humans = zeros(SVector{dwells[i].size,Int64})
    end
end
export init_dwellings

## randomly assigns age, age_group, preterm < 1 and smoker 19+
function apply_charectristics(x::Human)
    x.agegroup = rand(agedist)
    x.age = rand(agebraks[x.agegroup])

    if x.agegroup == 1   #assign preterm if < 1 years old
        x.preterm = rand(predist)
    end

    if x.agegroup == 7  #assign smoker if 19+ years old
        x.smoker = rand(smokdist)
    end
end
export apply_charectristics

###### Assign charectristics to humans: initial populations
function init_population()
    @inbounds for i = 1:length(humans)
        apply_charectristics(humans[i]) #charectristics; age,agegroup,preterm,smoker
    end
end
export init_population

######## randomly distributes humans to privet dwellings
function apply_humandistribution(x::Human, y::Dwelling)
    humans_G7 = findall(x -> x.agegroup == 7, humans) #find 19+ yearsold
    humans_G6 = findall(x -> x.agegroup == 6, humans) #find 5 > 19 yearsold
    humans_G5 = findall(x -> x.agegroup == 5, humans) #find 4:5 yearsold
    humans_G4 = findall(x -> x.agegroup == 4, humans) #find 3:4 yearsold
    humans_G3 = findall(x -> x.agegroup == 3, humans) #find 2:3 yearsold
    humans_G2 = findall(x -> x.agegroup == 2, humans) #find 1:2 yearsold
    humans_G1 = findall(x -> x.agegroup == 1, humans) #find < 1 yerasold

     for i in humans_G7
        #eligible dwells are not occuppied fully and at least has one free space for humans
        eligible_dwells = findall(y -> length(y.humans[y.humans .== 0]) > 0, dwells)
        id = rand(eligible_dwells)
        idd = findfirst(z -> z == 0, dwells[id].humans)
        dwells[id].humans[idd] = humans_G7[i]
    end
end
export apply_humandistribution



end # module
