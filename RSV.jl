## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh & Affan

module RSV

## Loading necessary packages
#using DataFrames

#using ClusterManagers
using Distributed
using Distributions
using StatsBase
using StaticArrays
using Random
using Parameters
using DataFrames
using CSV

#---------------------------------------------------
# Construct agent's types as an object
#---------------------------------------------------
mutable struct Human
    idx::Int64
    health::Int64       # 0 = susc, 1 = infected
    age::Int64          # in months
    agegroup::Int64    # G1 = 0-2, G2= 3-5, G3= 6-11, G4= 12-23, G5= 24-35, G6= 5:19 years, G7=19+ years
    preterm:: Bool     # true/false
    houseid::Int64
    Human() = new()
end

mutable struct Dwelling
    idx::Int64
    size::Int64      # household size
    avaisize::Int64  # left spots for humans
    infected::Int64  # infected=1 or not=0
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
const numInterestedAgeGroup = 5

### charectristics parameters
# categorical probability distribution of discrete age_groups
const agedist_Nuk =  Categorical(@SVector[0.006819215,0.007576906,0.008183058,0.026140324,0.026140324,0.026140324,0.022730717,0.310653129,0.565616002])
const agebraks_Nuk = @SVector[0:2, 3:5, 6:11, 12:23, 24:35, 36:47, 48:59, 60:227, 228:1200] #age_groups in months
const predist_Nuk = @SVector[0.033333333,0.09,0.021296296]  # distribution of preterm infants (0-2,3-5,6-11 months) obtained from data in file 'parameters_Shokoofeh.xlsx'
#const smokdist_Nuk = Categorical(@SVector[0.37, 0.63])  # distribution of adult smoker

## Dwelld Household Parameters
# total dwells per household size (0 is added for technical trick for writing a signle loop later)
const dwells_householdsize_Nuk = @SVector[0,755,595,555,585,1135] # householdsize of 1,2,3,4,5:20

# infect_dwells Parameters
const inf_range = @SVector[100:400]
#const SAR_Nuk = @SVector[0.0,0.0,0.0,0.0,1.0] # Nunavik: secondary attack rate for < 1 , 2, 3 yearsold kids
const SAR_Nuk = @SVector[0.63,0.63,0.63,0.40,0.27] # Nunavik: secondary attack rate for < 1 , 2, 3 yearsold kids
#Extreme Test: const SAR_Nuk = @SVector[1,1,1]

# Agents Arrays
const humans = Array{Human}(undef, population_Nuk)
const dwells = Array{Dwelling}(undef, dwells_Nuk)
export humans, dwells





#-----------------------------------------
#               SIMULATION FUNCTIONS
#-----------------------------------------

# run simulations for number of si and region "Nunavik or Nunavit"
function sim(si::Int64,reg::String)
    # collecting infection for different agegroups in df
    df = DataFrame(introduced_inf=Int64[], inf0to2=Int64[], inf3to5=Int64[], inf6to11=Int64[], inf12to23=Int64[], inf24to35=Int64[])
    @inbounds for x=1:si
        data = main(x,reg)
        push!(df,data)
    end
    return df
end
export sim

# producing average of infection for each agegroup
function ave_inf(df)
    nrows, ncols = size(df) #size of DataFrame
    # collecting average infection in ave
    ave = DataFrame(Name="Average", total_inf=0.0,inf0to2=0.0,inf3to5=0.0,inf6to11=0.0,inf12to23=0.0,inf24to35=0.0) # collecting data
    @inbounds for col =1:ncols
        ave[1,col+1] = sum(df[:,col])/nrows # average infections for each agegroup (among number of simulations)
    end
    return ave
end
export ave_inf

# saving raw results and averages
function save(df,ave)
    CSV.write("RSV_result_500sim.csv",df)
    CSV.write("RSV_average_500sim.csv",ave)
end
export save



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
    n_int = numInterestedAgeGroup

    ## simulation setup functions
    if simnumber == 1
        init_humans()
        init_dwellings(householdsize)
    else
        reset_humans(humans)
        reset_dwells(dwells)
    end
    init_population(humans,agedist,agebraks,predist)
    apply_humandistribution(dwells, humans)

    #infect dwells
    n_inf = num_inf(range)
    infect_dwells(dwells, humans, n_inf, SAR, n_int)

    ## colecting total infections for kids between 0-2,3-5,6-11,12,23,24-35 months
    ## data collection arrays
    infection_total = zeros(Int64, n_int+1) # total infected kids between 1-3 years old during one season
    infection_total[1] =n_inf
    @inbounds for i=1:n_int+1
        infection_total[i+1] = length(findall(x -> x.agegroup == i && x.health == 1 ,humans))
    end

    """
    ## run test for SAR =1 for only one specific agegroup
    ttt = findall(y-> y.==1, SAR)
    if length(ttt) > 0
        t = test(ttt[1])
    end
    return infection_total, t
    """

    # return number of total introduced infection, total infection occured for agegroups, run test
    return infection_total
end
export main

#--------------------------------------------------
# Required Functions for model
#---------------------------------------------------

## setting up agent 'humans'
function init_humans()
    @inbounds for i = 1:length(humans)
        humans[i] = Human()              ## create an empty human
        humans[i].idx = i
        humans[i].health = 0
        humans[i].age = 0
        humans[i].agegroup = 0
        humans[i].preterm = false
        humans[i].houseid = 0
    end
end
export init_humans

## randomly assigns age, age_group, preterm < 1
function apply_charectristics(x,agedist,agebraks,predist)
    _agegp = rand(agedist)::Int64
    x.agegroup = _agegp
    x.age = rand(agebraks[_agegp])

    #assign preterm if < 1 years old
    if _agegp == 1
        if rand() < predist[1]
            x.preterm = true
        end

    elseif _agegp == 2
        if rand() < predist[2]
            x.preterm = true
        end

    elseif _agegp == 3
        if rand() < predist[3]
            x.preterm = true
        end

    end
end
export apply_charectristics

###### Assign charectristics to humans: initial populations
function init_population(h,agedist,agebraks,predist)
    @inbounds for (i, element) in enumerate(humans)
        apply_charectristics(element,agedist,agebraks,predist)
    end
end
export init_population




## setup dwells with assigned humans and household size
function init_dwellings(dwells_householdsize)
    for i = 1:length(dwells)
        dwells[i] = Dwelling()  ## create an empty house
        dwells[i].idx = i
        dwells[i].infected = 0
    end
    S_min = 1
    S_max = 0
    sum_population = 0
    @inbounds for H = 1:length(dwells_householdsize)-1 # groupsize loop
        S_min += dwells_householdsize[H]
        S_max += dwells_householdsize[H+1]
        sum_population += H*dwells_householdsize[H+1]   ## total population up to the housesize = 5
        for S = S_min:S_max  # assigning size loop
            dwells[S].size = H
            dwells[S].avaisize = H
        end

        #### For last group size, randomly choose a size from household 5:20 (we dont have exact data for size 5>)
        if H == length(dwells_householdsize)-1
            popG9 = length(humans)-sum_population # population left to be distributaed for size 5>
            while popG9 !==0
                # shouldnt choose a size more than leftover people
                if popG9 >= 15
                    sizeG9 = rand(0:15)
                else
                    sizeG9 = rand(0:popG9)
                end
                G9 = rand(S_min:S_max) # dwells index must be between S_min and S_max
                dwells[G9].size += sizeG9
                dwells[G9].avaisize += sizeG9
                popG9 -= sizeG9
            end
        end
    end
end
export init_dwellings

###########
function find_a_house(d)
    available_house = findall(x -> x.avaisize >= 1, d)
    return rand(available_house)
end
export find_a_house

######## randomly distributes humans to privet dwellings
function apply_humandistribution(d, h)
    # first distribute one adult 19+ to every dwell
    humans_G9 = filter(x -> x.agegroup == 9 ,h)
    @inbounds for i = 1:length(d)
        humans_G9[i].houseid = i
        d[i].avaisize -= 1
    end

    # second distribute the leftover population, randomly into the dwells
    left_humans = filter(e -> !(e in humans_G9[1:length(d)]),h) #left humans to distribute
    @inbounds for j = 1:length(left_humans)
        idd = find_a_house(d)
        left_humans[j].houseid = idd
        d[idd].avaisize -= 1
    end
end
export apply_humandistribution


##================================================
# Introducing infection to the households
#-------------------------------------------------
## start infecting individuals who bring infection to the home
function num_inf(r)
    n_inf = rand(r[1])
    return n_inf
# Extreme Test: n_inf = length(dwells)
end
export num_inf

function infect_dwells(d, h, n_inf, SAR, n_int)
    ## input n_inf number of infected individuals between 5-100 yearsold. who bring infection homes
    dwell_inf = sample(d, n_inf; replace = false)  # infected dwells
    sec_att_rate = SAR
    @inbounds for i = 1:n_inf
        dwell_inf[i].infected = 1
        for j = 1:n_int
            sus = filter(x-> x.houseid == dwell_inf[i].idx && x.agegroup == j ,h) # susceptiple infants
            if length(sus) > 0
                for k = 1:length(sus)
                    if rand() < sec_att_rate[j]
                        sus[k].health = 1 # infected
                    end
                end
            end
        end
    end
end
export infect_dwells

#---------------------------------
# Reset Functions
#---------------------------------
# reset health and houseid for the next simulation
function reset_humans(h)
    @inbounds for i=1:length(h)
        h[i].health = 0
        h[i].age = 0
        h[i].agegroup = 0
        h[i].preterm = false
        h[i].houseid = 0
    end
end
export reset_humans

function reset_dwells(d)
    @inbounds for i=1:length(d)
        si = d[i].size
        d[i].avaisize = si
        d[i].infected = 0
    end
end
export reset_dwells

#------------------------------------------------------
# Test functions to see how good is model and the code
#-----------------------------------------------------

# make sure total population and dwells' sizes are equall
function population_equall_to_dwellSize(d,h)
    sum_size = 0
    sum_avaisize = 0
    l_h = length(h)
    @inbounds for (i,elem) in enumerate(d)
        sum_size += elem.size
        sum_avaisize += elem.avaisize
    end
    if sum_size !== l_h
        error("total population is NOT equall to total dwells' size")
    end
    if sum_avaisize !== 0
        error("There is at least one dwell with availible spot for more humans (Wrong)")
    end
    return println("total population is equall to total dwells' size"), println("Capacity of dwells are full")
end
export population_equall_to_dwellSize


## making sure all dwells include at least one adult (19+)
function at_least_one_adult(d)
    @inbounds for i = 1:length(d)
        dwells9 = filter(x -> x.houseid == i && x.agegroup == 9 , humans)
        if length(dwells9) == 0
            error("Not all dwells include at least one adult 19+")
        end
    end
    return println("All dwells include at least one adult 19+")
end
export at_least_one_adult



# see how many dwells with a specific agegroup=1or2or3... is infected?
function test(age)
    @inbounds a = findall(x -> x.infected == 1, dwells)
    tt = 0
    for i=1:length(a)
        if length(findall(x-> x.houseid == a[i] && x.agegroup == age, humans)) > 0
            tt += 1
        end
    end
    return tt
end
export test




#--------------------------------------
# to run in REPL
@everywhere using .RSV
####################
end # module

#using Distributed
#using ClusterManagers
#addprocs(SlurmManager(np);kwargs)
#@everywhere include ("RSV.jl")


#### test results ####
#----------------------------------------------------------------
# Row | infection | 0-2  | 3-5  | 6-11 | 12-23 | 24- 36 | test  |
#----------------------------------------------------------------
#  1  |    371    |  13  |  0   |  0   |   0   |   0    |  12   |
#  2  |    371    |  0   |  11  |  0   |   0   |   0    |  11   |
#  3  |    371    |  0   |  0   |  12  |   0   |   0    |  11   |
#  4  |    371    |  0   |  0   |  0   |   29  |   0    |  26   |
#  5  |    371    |  0   |  0   |  0   |   0   |  42    |  38   |
