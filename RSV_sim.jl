## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh

### Before run: include ("RSV_sim.jl") in Juno Editor

module RSV_sim

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
using Match


#---------------------------------------------------
# Construct agent's types as an object
#---------------------------------------------------
mutable struct Human
    idx::Int64
    health::Int64           # 0 = susc, 1 = infected
    age::Int64              # in months
    agegroup::Int64        # G1 = 0-2, G2= 3-5, G3= 6-11, G4= 12-23, G5= 24-35, G6= 5:19 years, G7=19+ years
    preterm_ill:: Bool     # true/false
    houseid::Int64
    Human() = new()
end
export Human

mutable struct Dwelling
    idx::Int64
    size::Int64      # household size
    avaisize::Int64  # left spots to include humans
    infected::Int64  # infected=1 or not=0
    Dwelling() = new()
end
export Dwelling

#----------------------------------------------------
# Global variables
#-----------------------------------------------------
# general parameters
# age populations are downloaded from Census Canada 2016
# data organized in file 'parameters_Shokoofeh.xlsx'
const numhumans = 13284 # total population of Nunavik
const numdwells = 3625 # total privet dwells in Nunavik
const n_int = 5

### charectristics parameters
# categorical probability distribution of discrete age_groups
const agedist =  Categorical(@SVector[0.006925625,0.007452575,0.014528756,0.025971093,0.025971093,0.025971093,0.022583559,0.308641975,0.561954231])
const agebraks = @SVector[0:2, 3:5, 6:11, 12:23, 24:35, 36:47, 48:59, 60:227, 228:1200] #age_groups in months
const predist = @SVector[0.093557786,0.076953948,0.046704722,  0.020289855]  # distribution of preterm_ill infants (0-2,3-5,6-11,12-23 months) obtained from data in file 'parameters_Shokoofeh.xlsx'



## Dwells Household Parameters
# total dwells per household size (0 is added for technical trick for writing a signle loop later)
const householdsize = @SVector[0,755,595,555,585,1135] # householdsize of 1,2,3,4,5:20

# Nunavik: mean secondary attack rate for < 1 , 2, 3 yearsold kids from Kenya study:
#Patrick K. Munywoki et. al., The Source of Respiratory Syncytial Virus Infection In Infants:
#A Household Cohort Study In Rural Kenya, The Journal of Infectious Diseases,
# https://doi.org/10.1093/infdis/jit828
#const SAR = @SVector[0.0,0.0,0.0,0.0,1.0]
const SAR = @SVector[0.636,0.636,0.636,0.40,0.27]
#Extreme Test: const SAR = @SVector[1,1,1]
const ns = 44 # number of sampled households in Kenya study

# Agents Arrays
const humans = Array{Human}(undef,numhumans)
const dwells = Array{Dwelling}(undef, numdwells)
export humans, dwells


#--------------------------------------------------------------
#               SIMULATION FUNCTION
#--------------------------------------------------------------

# run simulations for number of si and intensity of season
#sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
# example to run: sim(500,:mild)
function sim(si::Int64, season::Symbol)
    # collecting infection for different agegroups in df
    names_inf = Symbol.(["introduced_inf","inf0to2_total","inf0to2_H","inf0to2_P",
     "inf3to5_total","inf3to5_H","inf3to5_P","inf6to11_total","inf6to11_H","inf6to11_P",
      "inf12to23_total","inf12to23_H","inf12to23_P","inf24to35"])
    df_infection = DataFrame([Int64 for i=1:length(names_inf)], names_inf)

    #collecting population per agegroup
    names_pop = Symbol.(["pop0to2_total","pop0to2_H","inf0to2_P",
     "pop3to5_total","pop3to5_H","pop3to5_P","pop6to11_total","pop6to11_H","pop6to11_P",
      "pop12to23_total","pop12to23_H","pop12to23_P","pop24to35"])
    df_population = DataFrame([Int64 for i=1:length(names_pop)], names_pop)

    # collecting incidence rates per agegroups and conditions, i.e. healthy and preterm_ill
    names_inc = Symbol.(["introduced_inf","inc0to2_H","inc0to2_P","inc3to5_H","inc3to5_P",
    "inc6to11_H","inc6to11_P","inc12to23_H","inc12to23_P","inc24to35"])
    df_incidence = DataFrame([Float64 for i=1:length(names_inc)], names_inc)

    n=zeros(si) #to collect number of dwells with one infant < 12 months
    # collecting incidence rates for regional hospital admission



    @inbounds for x=1:si #simulation
        inf_data, pop_data, inc_data = main(x,season) # calculate infection and population per simulation
        # include data into dataframe for saving format
        push!(df_infection,inf_data)
        push!(df_incidence,inc_data)
        if season == :mild
            push!(df_population,pop_data)
        end
    end

    #---------
    # saving raw results for 500 simulations
    CSV.write("infection_$season.csv",df_infection)
    CSV.write("incidence_$season.csv",df_incidence)
    if season == :mild
        CSV.write("population.csv",df_population)
    end
    return println("$season season is done")
end
export sim




#-------------------------------------------------
# main function to run for one simulation
#-------------------------------------------------
function main(simnumber::Int64, season::Symbol)
    Random.seed!(simnumber)

    ## simulation setup functions
    if simnumber == 1
        init_humans()
        init_dwellings(householdsize)
    else
        reset_humans(humans)
        reset_dwells(dwells)
    end

    init_population(humans,agedist,agebraks,predist)
    #apply_dwellsSize(dwells,humans,householdsize)
    apply_humandistribution(dwells, humans)

    #--------------infecting dwells
    # find dwells with minimum one humans < 1 yearsold
    dwells_to_infect = find_dwells_to_infect()
    n = length(dwells_to_infect) # total number of eligible dwells to infect
    range = range_of_inf(n,season) # get range of outbreak depending on season

    # total number of households introduced to RSV per simulation
    n_inf = rand(range)
    infect_dwells(dwells, dwells_to_infect, humans, n_inf, SAR, ns, n_int)

    #-----------------
    #------ Results (objective agegroups are children under 3 years old)
    #----------------
    #caculate total infection and population per agegroups
    infection, population = calcu_infection_population(n_int)

    # Choose hosp incidence rates for every agegroup per simulation
    incidence = choose_incidence(n_int,season)

    # add initial infection to the first coloumn
    pushfirst!(infection, n_inf)
    pushfirst!(incidence, n_inf)




    """
    ## run test for SAR =1 for only one specific agegroup
    ttt = findall(y-> y.==1, SAR)
    if length(ttt) > 0
        t = test(ttt[1])
    end
    return infection_crt, t
    """
    return infection, population, incidence
end
export main

#--------------------------------------------------
# Required Functions for model
#---------------------------------------------------

## setting up agent 'humans'
function init_humans()
    @inbounds for i = 1:length(humans)
        humans[i] = Human()              ## create an empty human
        x = humans[i]
        x.idx = i
        x.health = 0
        x.age = 0
        x.agegroup = 0
        x.preterm_ill = false
        x.houseid = 0
    end
end
export init_humans

###### randomly assigns age, age_group, underlying condition to children < 2 years
## underlying conditions include prematurity born < 1 years and chronically-ill children < 2 years old
function apply_charectristics(x,agedist,agebraks,predist)
    _agegp = rand(agedist)::Int64
    x.agegroup = _agegp
    x.age = rand(agebraks[_agegp])

    #assign preterm_ill if < 1 years old
    if _agegp == 1
        if rand() < predist[1]
            x.preterm_ill = true
        end

    elseif _agegp == 2
        if rand() < predist[2]
            x.preterm_ill = true
        end

    elseif _agegp == 3
        if rand() < predist[3]
            x.preterm_ill = true
        end
    elseif _agegp == 4
        if rand() < predist[4]
            x.preterm_ill = true
        end

    end
end
export apply_charectristics

###### Assign charectristics to humans: initial populations
function init_population(h,agedist,agebraks,predist)
    @inbounds for (i, element) in enumerate(h)
        apply_charectristics(element,agedist,agebraks,predist)
    end
end
export init_population



"""
## setup houses called dwells with initial charectristics
function init_dwellings()
    @inbounds for i = 1:length(dwells)
        dwells[i] = Dwelling()  ## create an empty house
        dwells[i].idx = i
        dwells[i].size = 0      ## size of households
        dwells[i].avaisize = 0  ##
        dwells[i].infected = 0  ## infected or not?

    end
end
export init_dwellings

### Assign size to empty dwells
function apply_dwellsSize(d,h,householdsize)
    sum_population = 0
    num_size = length(householdsize) # number of household's distribution for size=1, 2, 3, 4, 5to20
    @inbounds for H = 1:num_size
        dwells_need_size = filter(x-> x.size == 0, d) # remaining dwells to be assigned size
        # randomly, sample dwells for householdsize H
        dwellsSize =sample(dwells_need_size,householdsize[H]; replace=false)
        for S=1:length(dwellsSize)
            dwellsSize[S].size = H #assign size H=1,2,3,4,5
            dwellsSize[S].avaisize = H #assign availiblity for humans
        end
        ## add-up population, while assigning size. We need to normalize total population to the dwells size
        sum_population += H*householdsize[H]

        # assigning size H=5to20
        if H == num_size
            pop5to20 = length(h)-sum_population # remained population to be normalized to household size 5:20
            while pop5to20 !==0
                if pop5to20 >= 15
                    Si = rand(0:15)
                else
                    Si = rand(0:pop5to20)
                end
                dwells5to20 = sample(dwellsSize)
                dwells5to20.size += Si #assign size H=1,2,3,4
                dwells5to20.avaisize += Si #assign availiblity for humans
                pop5to20 -= Si
            end
        end
    end
end
export apply_dwellsSize
"""


function init_dwellings(dwells_householdsize)
    for i = 1:length(dwells)
        dwells[i] = Dwelling()  ## create an empty house
        x = dwells[i]
        x.idx = i
        x.infected = 0
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
                #if dwells[G9].size >= 6 # shouldnt select previously chosen dwells (more than 5 means it was given size before)
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

function get_percentage_of_outbreak(sea::Symbol)
    ret = @match sea begin
        :mild => StaticArrays.@SVector[0.3,0.5]
        :moderate => StaticArrays.@SVector[0.5,0.7]
        :severe => StaticArrays.@SVector[0.7,0.9]
        #_ => error("what is the season?")
    end
    return ret
end
export get_percentage_of_outbreak


function range_of_inf(n::Int64,season::Symbol)
    #input is number of dwells with infants <1 yearsold
    p = get_percentage_of_outbreak(season) #mild/moderate/severe season?
    r1 = floor(Int,p[1]*n)
    r2 = floor(Int,p[2]*n)
    range = r1:r2
    # gives back number of eligiblie dwells to infect
    # and the 30 to 50 % range of dwell outbreaks
    return range
end
export range_of_inf



#-------------Find eligible dwells to infect
function find_dwells_to_infect()
    # find children under 1 years old
    children_under1years = filter(x-> x.agegroup < 4 ,humans)
    # find houses of which these infants < 1 yearold residing in it
    dwells_to_infect = [children_under1years[i].houseid for i=1:length(children_under1years)]
    unique(dwells_to_infect) # remove duplicate houses
    return dwells_to_infect
end
export find_dwells_to_infect

## -------INFECT DWELLS, which at least include one children under 3 yearsold
## definition of inputs:
# d:dwells, h:humans, n_inf: initial infection given to dwells,
# SAR:secondary attack rates, ns: number of samples for SAR
#n_int:number of objective agegroup we are interested in
#-----------------------------------------------------------
function infect_dwells(d, d_to_infect, h, n_inf, SAR, ns, n_int)
    dwells_inf = sample(d_to_infect, n_inf; replace = false)  # infected dwells
    @inbounds for i = 1:n_inf
        infected_id = dwells_inf[i]
        d[infected_id].infected = 1
        for j = 1:n_int
            sus = filter(x-> x.houseid == infected_id && x.agegroup == j ,h) # susceptiple infants
            if length(sus) > 0
                for k = 1:length(sus)
                    if rand() < rand(Beta(SAR[j]*ns, (1-SAR[j])*ns))
                        sus[k].health = 1 # infected
                    end
                end
            end
        end
    end
end
export infect_dwells


#--------------------------------------------------------------
#                   Infection & Population
#-------------------------------------------------------------
function calcu_infection_population(n_int::Int64)
    # data-collection vectors
    infection = zeros(Int64, n_int*3-2)
    population = zeros(Int64,n_int*3-2)

    # iteration for our collecting arrays
    itr =@SVector[1,4,7,10,13]
    @inbounds for i=1:n_int
        j = itr[i]
        # filter objective population in agegroup 'i'
        humans_agegroup_i = filter(x-> x.agegroup == i, humans)
        pop_t = length(humans_agegroup_i) #total population per agegroup

        if i !== 5 # kids under than 2 yearsold
            # population (total,healthy,preterm)
            pop_h = length(findall(y-> y.preterm_ill == false , humans_agegroup_i)) #healthy
            pop_p = pop_t - pop_h #preterm_ill

            # collect population data
            #-----------------------------------
            population[j] = pop_t #total pop per agegroup
            population[j+1] = pop_h # pop per healthy infants
            population[j+2] = pop_p # pop per preterm_ill infants

            # calculate total infection
            #-----------------------------------
            # infected kids (total,healthy,preterm)
            t = length(findall(z -> z.health == 1, humans_agegroup_i)) #total infected
            h = length(findall(m -> m.health == 1 && m.preterm_ill == false, humans_agegroup_i)) #healthy infected
            p = t - h # preterm infected
            # collecting infection data
            infection[j] = t #total infection per agegroup
            infection[j+1] = h # infections per healthy infants
            infection[j+2] = p #infections per preterm_ill infants
        else
            ## 3 yearsold (NO preterm/healthy condition applied)
            # collect population
            #-------------------
            population[j] = pop_t #total pop per 3 yearsold

            # collect infection
            #------------------
            _h = length(findall(a -> a.health == 1, humans_agegroup_i))
            infection[j] = _h # total infection per 3 yearsold
        end
    end
    return infection, population
end
export calcu_infection_population




#--------------------------------------------------------------------------
# Regional Hosp Incidence Rates (initial incidence)
#-------------------------------------------------------------------------
# randomly choose an incidence rate for regional hospital admissions in every simulation (500 run) from given data range (clinical information)
# input; n_int::number of agegroups, n_inf:: number of introduced infected dwells
function choose_incidence(n_int::Int64, season::Symbol)
    # data-collection vectors
    #inc = zeros(Float64,n_int*2-1) # hosp incidence rates
    inc = zeros(Float64,n_int*2-1)

    # get incidence range which depending on the season (mild/moderate/severe)
    rate = get_incidence_of_season(season)

    # iteration for our collecting arrays
    itr =@SVector[1,3,5,7,9]
    @inbounds for i in itr
        # selecting random incidence rate
        inc[i] = rand(rate[i]) #healthy
        inc[i+1] = rand(rate[i+1]) #preterm_ill
    end
    return inc
end
export choose_incidence


# Hospitalization rates per 100,000 person-time for Nunavik (provided by collaborator)
function get_incidence_of_season(season::Symbol)
    s = @match season begin
        #rates corrospond to our agegroups G1 [healthy,preterm_ill] to G5 ****(G4 and G5 are missing for now)?????
        # G1:h,p - G2:h,p - G3:h,p - G4:h,p - G5: total
        :mild => StaticArrays.@SVector[27.1:50.5, 30.69:139.34, 37.2:41.2, 213.68:265.23, 20.6:28.6, 207.58:254.28, 0.0:0.0, 0.0:0.0, 0.0:0.0]
        :moderate => StaticArrays.@SVector[50.6:74,139.35:248, 41.3:45.3, 265.24:316.79, 28.7:36.7, 254.29:300.99, 0.0:0.0, 0.0:0.0, 0.0:0.0]
        :severe => StaticArrays.@SVector[74.1:97.3,248.01:356.65, 45.4:49.4, 316.80:368.32, 36.8:44.8, 301:347.7, 0.0:0.0, 0.0:0.0, 0.0:0.0]
    end

    return s
end
export get_incidence_of_season

#---------------------------------
# Reset Functions
#---------------------------------
# reset health and houseid for the next simulation
function reset_humans(h)
    for i=1:length(h)
        x = h[i]
        x.health = 0
        x.age = 0
        x.agegroup = 0
        x.preterm_ill = false
        x.houseid = 0
    end
end
export reset_humans

function reset_dwells(d)
    for i=1:length(d)
        x = d[i]
        si = x.size
        x.avaisize = si
        x.infected = 0
    end
end
export reset_dwells


#---------------------------------------------
# producing average of df_data for each coloumn
function get_average(df)
    nrows, ncols = size(df) #size of DataFrame
    # collecting average infection in ave
    ave = zeros(Int64, ncols)
    @inbounds for col =1:ncols
        #ave[1,col+1] = sum(df[:,col])/nrows # average (among number of simulations)
        ave[col] = round(sum(df[:,col])/nrows) # average (among number of simulations)
    end
    return ave
end
export get_average



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
@everywhere using .RSV_sim
####################
end # module

#using Distributed
#using ClusterManagers
#addprocs(SlurmManager(np);kwargs)
#@everywhere include ("RSV_sim.jl")


#### test results ####
#----------------------------------------------------------------
# Row | infection | 0-2  | 3-5  | 6-11 | 12-23 | 24- 36 | test  |
#----------------------------------------------------------------
#  1  |    371    |  13  |  0   |  0   |   0   |   0    |  12   |
#  2  |    371    |  0   |  11  |  0   |   0   |   0    |  11   |
#  3  |    371    |  0   |  0   |  12  |   0   |   0    |  11   |
#  4  |    371    |  0   |  0   |  0   |   29  |   0    |  26   |
#  5  |    371    |  0   |  0   |  0   |   0   |  42    |  38   |
