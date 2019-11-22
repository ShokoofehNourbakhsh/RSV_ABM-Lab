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
    agegroup::Int64    # G1 = < 1, G2= 1, G3= 2, G4= 3, G5= 4, G6= 5:19, G7= > 19
    smoker::Int64       # 1 = no,  2 = yes
    preterm:: Int64     # 1 = no,  2 = yes
    Human() = new()
end

mutable struct Dwelling
    idx::Int64
    humans::Array{Int64}
    #humans::Human
    size::Int64
    Dwelling() = new()
end

######## global variables
# Nunavik age population. Downloaded from Census Canada 2016
# data organized in file 'parameters_Shokoofeh'
const agedist =  Categorical(@SVector[0.022727273, 0.026136364, 0.026136364, 0.026136364, 0.022727273, 0.310606061, 0.565530303]) # categorical probability distribution of discrete age_groups
const agebraks = @SVector[0:0.99, 1:1.99, 2:2.99, 3:3.99, 4:4.99, 5:18.99, 19:100] #age_groups
const predist = Categorical(@SVector[0.972027972, 0.027972028])  # distribution of preterm infants
const smokdist = Categorical(@SVector[0.37, 0.63])  # distribution of adult smoker
#const P = ModelParameters() # I'll make the parameter file
const population_Nuk = 13200 # total population of Nunavik
const dwells_Nuk = 3625 # total privet dwells in Nunavik
## household size categorizes into 1 person , two persons ... more than 5 persons in each household
const householdsize = @SVector[1,2,3,4,5:20]
const dwells_householdsize = @SVector[0,755,595,555,585,1135] # total dwells per household size (0 is added for technical trick for writing a signle loop later)
const humans = Array{Human}(undef, population_Nuk )
const dwells = Array{Dwelling}(undef, dwells_Nuk)


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

"""
######### setup dwells with assigned humans and household size
function init_dwellings()
    @inbounds for i = 1:length(dwells)
        dwells[i] = Dwelling()  ## create an empty house
        dwells[i].idx = i
    end
    S_min = 1
    S_max = 0
    for H = 1:length(dwells_householdsize)-1 # groupsize loop
        S_min += dwells_householdsize[H]
        S_max += dwells_householdsize[H+1]
        for S = S_min:S_max  # assigning size loop
            dwells[S].size = H
            dwells[S].humans = zeros(SVector{dwells[S].size,Int64})
        end
    end
    test_dwellsize = sum([dwells[n].size for n=1:length(dwells)]) # It must give total humnas
end
export init_dwellings, test_dwellsize

"""


######### setup dwells with assigned humans and household size
function init_dwellings()
    @inbounds for i = 1:length(dwells)
        dwells[i] = Dwelling()  ## create an empty house
        dwells[i].idx = i
    end
    S_min = 1
    S_max = 0
    sum_population = 0
    for H = 1:length(dwells_householdsize)-1 # groupsize loop
        S_min += dwells_householdsize[H]
        S_max += dwells_householdsize[H+1]
        sum_population += H*dwells_householdsize[H+1]   ## total population up to the housesize = 5
        for S = S_min:S_max  # assigning size loop
            dwells[S].size = H
            dwells[S].humans = zeros(SVector{dwells[S].size,Int64})
        end

        #### For last group size, randomly choose a size from household 5:20 (we dont have exact data for size 5>)
        if H == length(dwells_householdsize)-1
            popG5 = length(humans)-sum_population # population left to be distributaed for size 5>
            while popG5 !==0
                # shouldnt choose a size more than leftover people
                if popG5 >= 15
                    sizeG5 = rand(0:15)
                else
                    sizeG5 = rand(0:popG5)
                end
                G5 = rand(S_min:S_max) # dwells index must be between S_min and S_max
                dwells[G5].size += sizeG5
                dwells[G5].humans = zeros(SVector{dwells[G5].size,Int64})
                popG5 -= sizeG5
            end
        end
    end

    return
    if sum([dwells[n].size for n=1:length(dwells)]) == length(humans)
        println("total population is equall to total dwells' size")
    else
        println("total population is NOT equall to total dwells' size")
    end
end
export init_dwellings


######## randomly distributes humans to privet dwellings
function apply_humandistribution()
    # first distribute one adult 19+ to every dwell
    humans_G7 = findall(x -> x.agegroup == 7, humans)
    for i = 1:length(dwells)
        dwells[i].humans[1] = humans_G7[i]
    end
    # second distribute the leftover G7 and all other ages, randomly
    for j = 1:length(agebraks)
        if j == length(agebraks) # agegroup 7 (19+)
            humans_G = humans_G7[length(dwells)+1:end] # leftover G7
        else
            humans_G = findall(x -> x.agegroup == j, humans) # finding agegroups j=1,2,3,4,5:18,19>
        end
        for i = 1:length(humans_G)
            #eligible dwells are not occuppied fully and at least has one free space for humans
            eligible_dwells = findall(y -> length(y.humans[y.humans .== 0]) > 0, dwells)
            if  length(eligible_dwells) == 0 #  should never happen, otherwide our dwells.humans not setup correctly
                error("no house left to distribute humans")
            end
            id = rand(eligible_dwells)
            idd = findfirst(z -> z == 0, dwells[id].humans)
            dwells[id].humans[idd] = humans_G[i]
        end
    end
    ## making sure dwells.size is setup correctly
    eligible_dwells = findall(y -> length(y.humans[y.humans .== 0]) > 0, dwells)
    if  length(eligible_dwells) !== 0 #
        error("total househols's size is larger than total population")
    end
    ## making sure all dwells include at least one adult (19+)
    for k = 1:length(dwells)
        dwells7 = findall(t -> humans[t].agegroup == 7 , dwells[k].humans)
        if length(dwells7) == 0
            error("Not all dwells include at least one adult 19+")
        end
    end
end
export apply_humandistribution

##### introducing infection to the households
function infect_dwells(n_inf::Int64)
    ## input n_inf number of infected individuals between 5-100 yearsold. who bring infection home
    dwell_inf = sample(dwells, n_inf; replace = false)  # infected dwells
    for i = 1:n_inf
        sus_1 = findall(x -> humans[x].agegroup == 1 ,dwell_inf[i].humans) # susceptiple infants
        sus_2 = findall(x -> humans[x].agegroup == 2 ,dwell_inf[i].humans) # susceptiple 1-2 years
        sus_3 = findall(x -> humans[x].agegroup == 3 ,dwell_inf[i].humans) # susceptiple 2-3 years
        if length(sus_1) == 0 && length(sus_2) == 0 && length(sus_3) == 0
            continue
        elseif length(sus_1) > 0
            for j1 = 1:length(sus_1)
                r = rand(Categorical(@SVector[0.38,0.62]))
                if r == 1
                    humans[sus_1[j1]].health = 0  # healthy
                else
                    humans[sus_1[j1]].health = 1 # infected
                end
            end
        elseif length(sus_2) > 0
            for j2 = 1:length(sus_2)
                r = rand(Categorical(@SVector[0.6,0.4]))
                if r == 1
                    humans[sus_2[j2]].health = 0  # healthy
                else
                    humans[sus_2[j2]].health = 1 # infected
                end
            end
        elseif length(sus_3) > 0
            for j3 = 1:length(sus_3)
                r = rand(Categorical(@SVector[0.73,0.27]))
                if r == 1
                    humans[sus_3[j3]].health = 0  # healthy
                else
                    humans[sus_3[j3]].health = 1 # infected
                end
            end
        end
    end
end
export infect_dwells



end # module
