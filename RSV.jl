## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh & Affan

module RSV

## Loading necessary packages
#using DataFrames
using Distributions
using StatsBase
using StaticArrays
using Random
using Parameters

#---------------------------------------------------
# Construct agent's types as an object
#---------------------------------------------------
mutable struct Human
    idx::Int64
    health::Int64       # 0 = susc, 1 = infected
    age::Int64          # in months
    agegroup::Int64    # G1 = 0-2, G2= 3-5, G3= 6-11, G4= 12-23, G5= 24-35, G6= 5:19 years, G7=19+ years
    preterm:: Int64     # 1 = no,  2 = yes
    Human() = new()
end

mutable struct Dwelling
    idx::Int64
    humans::Array{Int64}
    size::Int64
    Dwelling() = new()
end



#----------------------------------------------------
# Global variables
#-----------------------------------------------------
# general parameters
# age populations are downloaded from Census Canada 2016
# data organized in file 'parameters_Shokoofeh.xlsx'
const population_Nuk = 13198 # total population of Nunavik
const dwells_Nuk = 3625 # total privet dwells in Nunavik
const population_Nut = 1500 #???????
const dwells_Nut = 5000 #???????
const region = ["Nunavik","Nunavit"]

### charectristics parameters
# categorical probability distribution of discrete age_groups
const agedist_Nuk =  Categorical(@SVector[0.006819215,0.007576906,0.008183058,0.026140324,0.026140324,0.026140324,0.022730717,0.310653129,0.565616002])
const agebraks_Nuk = @SVector[0:2, 3:5, 6:11, 12:23, 24:35, 36:47, 48:59, 60:227, 228:1200] #age_groups in months
#const agebraks_Nuk = @SVector[0:0.99, 1:1.99, 2:2.99, 3:3.99, 4:4.99, 5:18.99, 19:100] #age_groups
#const agedist_Nuk =  Categorical(@SVector[0.022727273, 0.026136364, 0.026136364, 0.026136364, 0.022727273, 0.310606061, 0.565530303])
const predist_Nuk = @SVector[0.033333333,0.09,0.021296296]  # distribution of preterm infants (0-2,3-5,6-11 months)
#const smokdist_Nuk = Categorical(@SVector[0.37, 0.63])  # distribution of adult smoker

## dwell household parameters
const householdsize_Nuk = @SVector[1,2,3,4,5:20]
# total dwells per household size (0 is added for technical trick for writing a signle loop later)
const dwells_householdsize_Nuk = @SVector[0,755,595,555,585,1135]

# infect_dwells parameters
const inf_range = @SVector[100:400]
const SAR_Nuk = @SVector[0.63,0.63,0.63,0.40,0.27] # Nunavik: secondary attack rate for < 1 , 2, 3 yearsold kids
#Extreme Test: const SAR_Nuk = @SVector[1,1,1]

const humans = Array{Human}(undef, population_Nuk)
const dwells = Array{Dwelling}(undef, dwells_Nuk)
export humans, dwells, infection_total



#-------------------------------------------------
# main function to run for one simulation
#-------------------------------------------------
function main(simnumber::Int64 , region::String)
    Random.seed!(simnumber)

    # set the global parameters of the model
    if region == "Nunavik"
        Reg = region[1]
        numhumans  = population_Nuk
        agedist = agedist_Nuk
        agebraks = agebraks_Nuk
        predist = predist_Nuk
        numdwells = dwells_Nuk
        householdsize = dwells_householdsize_Nuk
        SAR = SAR_Nuk
    elseif region == "Nunavit"
        Reg = region[2]
        numhumans  = population_Nut
        numdwells = dwells_Nut
        SAR = SAR_Nut
    end
    range = inf_range



    ## data collection arrays
    infection_total = zeros(Int64, 5) # total infected kids between 1-3 years old during one season

    ## simulation setup functions
    if simnumber == 1
        init_humans(numhumans)
        init_dwellings(numdwells, householdsize)
    else
        reset_health(humans)
        reset_humandistribution(dwells)
    end
    init_population(humans, agedist, agebraks, predist)
    apply_humandistribution(dwells, humans)

    ## start infecting individuals who bring infection to the home
    n_inf = rand(range[1])
    # Extreme Test: n_inf = length(dwells)
    infect_dwells(dwells, humans, n_inf, SAR)

    ## colecting total infections for kids between 0-2,3-5,6-11,12,23,24-35 months
    for i=1:5
        infection_total[i] = length(findall(x -> x.agegroup == i && x.health == 1 ,humans))
    end
    return n_inf, infection_total
end
export main


#--------------------------------------------------
# Required Functions for model
#---------------------------------------------------

## setting up agent 'humans'
function init_humans(num_humans::Int64)
    @inbounds for i = 1:length(humans)
        humans[i] = Human()              ## create an empty human
        humans[i].idx = i
        humans[i].health = 0
        humans[i].age = 0
        humans[i].agegroup = 0
        humans[i].preterm = 1
    end
end
export init_humans

## randomly assigns age, age_group, preterm < 1
function apply_charectristics(x, agedist , agebraks, predist)
    x.agegroup = rand(agedist)
    x.age = rand(agebraks[x.agegroup])

    #assign preterm if < 1 years old
    if x.agegroup == 1
        x.preterm = rand(Categorical(@SVector[1-predist[1],predist[1]]))

    elseif x.agegroup == 2
        x.preterm = rand(Categorical(@SVector[1-predist[2],predist[2]]))

    elseif x.agegroup == 3
        x.preterm = rand(Categorical(@SVector[1-predist[3],predist[3]]))
    end
end
export apply_charectristics

###### Assign charectristics to humans: initial populations
function init_population(h, agedist , agebraks, predist)
    @inbounds for i = 1:length(h)
        apply_charectristics(h[i], agedist, agebraks, predist) #charectristics; age,agegroup,preterm
    end
end
export init_population




## setup dwells with assigned humans and household size
function init_dwellings(num_dwells::Int64, dwells_householdsize)
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
    if sum([dwells[n].size for n=1:length(dwells)]) == length(humans)
        a = println("total population is equall to total dwells' size")
    else
        a = println("total population is NOT equall to total dwells' size")
    end
    return a
end
export init_dwellings


######## randomly distributes humans to privet dwellings
function apply_humandistribution(d, h)
    # first distribute one adult 19+ to every dwell
    humans_G9 = filter(x -> x.agegroup == 9 ,h)
    @inbounds for i = 1:length(d)
        d[i].humans[1] = humans_G9[i].idx
    end

    # second distribute the leftover population, randomly into the dwells
    left_humans = filter(e -> !(e in humans_G9[1:length(d)]),h) #left humans to distribute
    for j = 1:length(left_humans)
        #eligible dwells are not occuppied fully and at least has one free space for humans
        eligible_dwells = findall(y -> length(y.humans[y.humans .== 0]) > 0, d)
        if  length(eligible_dwells) == 0 #  should never happen, otherwide our dwells.humans not setup correctly
            error("no house left to distribute humans")
        end
        id = rand(eligible_dwells)
        idd = findfirst(z -> z == 0, d[id].humans)
        d[id].humans[idd] = left_humans[j].idx
    end

    ## making sure dwells.size is setup correctly
    eligible_dwells = findall(y -> length(y.humans[y.humans .== 0]) > 0, d)
    if  length(eligible_dwells) !== 0 #
        error("total househols's size is larger than total population")
    end
    ## making sure all dwells include at least one adult (19+)
    for k = 1:length(d)
        dwells9 = findall(t -> h[t].agegroup == 9 , d[k].humans)
        if length(dwells9) == 0
            error("Not all dwells include at least one adult 19+")
        end
    end
end
export apply_humandistribution


##================================================
# Introducing infection to the households
#-------------------------------------------------
function infect_dwells(d, h, n_inf, SAR)
    ## input n_inf number of infected individuals between 5-100 yearsold. who bring infection homes
    dwell_inf = sample(d, n_inf; replace = false)  # infected dwells
    for i = 1:n_inf
        for j = 1:5
            sus = filter(x -> h[x].agegroup == j ,dwell_inf[i].humans) # susceptiple infants
            if length(sus) > 0
                for k = 1:length(sus)
                    if rand() < SAR[j]
                        h[sus[k]].health = 1 # infected
                    end
                end
            end
        end
    end
end
export infect_dwells


#---------------------------------
# Reset Funtions
#---------------------------------
function reset_health(h)
    ### reset health's status of kids for the next simulation
    kids = filter(x -> x.agegroup == 1 || x.agegroup == 2 || x.agegroup == 3 || x.agegroup == 4 || x.agegroup == 5,h)
    @inbounds for i=1:length(kids)
        kids[i].health = 0
    end
end
export reset_health

function reset_humandistribution(d)
    @inbounds for i=1:length(d)
        for j=1:d[i].size
            d[i].humans[j] = 0
        end
    end
end
export reset_humandistribution


####################
end # module
