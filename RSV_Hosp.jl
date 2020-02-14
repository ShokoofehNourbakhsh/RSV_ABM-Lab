## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh & Affan

module RSV_Hosp

## Loading necessary packages

#using ClusterManagers
using Distributed
using Distributions
#using StatsBase
using StaticArrays
#using Parameters
using DataFrames
using CSV
using Statistics
using Random


#----------------------------------------------------
# Global variables
#-----------------------------------------------------
# general parameters
const population_Nuk = 13216 # total population of Nunavik
const dwells_Nuk = 3625 # total privet dwells in Nunavik





#----------------------------------------------------------------------------
# Estimate Hospital Admissions and clinc outpatients for All senarios; 1,2,3,4,5
#-----------------------------------------------------------------------------
# S represent senarios;
# 1- No intervention 2- Palivizumab program < 2017, 3- Palivizumab program > 2017
# 4- Long-acting monoclonal antibody (LAMA), 5- Maternal vaccine
function HospClinic(S::Int64)

    #---------------------------
    # Reading total infection (results) from Agent-based Simulation
    # introduced_inf | inf0to2_total | inf0to2_H | inf0to2_P | ...500×14 DataFrame table
    infection = read_file("infection.csv") # infections per simulation in a row
    _inf = select(infection, Not([1,2,5,8,11])) # take out 'total infection' columns per agegroup
    inf = Matrix(_inf) # convert to matrix

    #--------------------------
    # senario 1; No intervention
    if S == 1
        ###### Reading Data Results from Agent-based Simulation
        # pop0to2_total | pop0to2_H | pop0to2_P | ... 500×13 DataFrame table
        population = read_file("population.csv") # population per simulation in  a row
        _pop = select(population, Not([1,4,7,10])) #remove columns (we only need healthy and preterm_ill population)
        pop = Matrix(_pop)

        # introduced_inf | inc0to2_H | inc0to2_P | ...  500X10 DataFrame table
        incidence =  read_file("incidence.csv") # hosp. incidence rates per agegroup
        _inc = select(incidence, Not(1)) #take out initial infection column
        inc = Matrix(_inc)

        #-----------------------------
        # calculate hospitalization rate's denominator per simulation
        # denom = person*days/100000
        denom = calcu_denominator(pop)

        #------------------------------
        # Determine hospitalization cases and outpatient numbers per every simulation
        # NO INTERVENTION (senario 1:S1)
        Hospcases, outpatients = calcu_HospClinicS1(inf,denom,inc)


    #--------------------------------------------
    #Senario 2; prophylaxis Palivizumab before 2017
    elseif S == 2

        #--------------------
        # Reading Data Results from HospS1 & clinicS1 (NO Intervention)
        # HospS1_0to2_H | HospS1_0to2_p |... 500×9 DataFrame table
        _HospS1 = read_file("HospS1.csv") # Hosp cases per simulation
        HospS1 = Matrix(_HospS1)

        #---------------------
        # Determine hospitalization cases and outpatient numbers per every simulation
        # S2: Palivizumab for preterm_ill, only. (senario 2)
        Hospcases, outpatients = calcu_HospClinicS2(inf,HospS1)
    end

    #----------------------------------------------------------
    # DataFrame to collect data
    # names
    names_cases = Symbol.(["HospS$S 0to2_H", "HospS$S 0to2_P","HospS$S 3to5_H", "HospS$S 3to5_P",
    "HospS$S 6to11_H","HospS$S 6to11_P","HospS$S 12to23_H", "HospS$S 12to23_P", "HospS$S 24to35"])

    names_patients = Symbol.(["clinicS$S 0to2_H", "clinicS$S 0to2_P","clinicS$S 3to5_H", "clinicS$S 3to5_P",
    "clinicS$S 6to11_H","clinicS$S 6to11_P","clinicS$S 12to23_H", "clinicS$S 12to23_P", "clinicS$S 24to35"])

    df_Hosp = DataFrame(Hospcases, names_cases)
    df_clinic = DataFrame(outpatients, names_patients)

    #---------
    # saving data
    CSV.write("HospS$S.csv",df_Hosp)
    CSV.write("clinicS$S.csv",df_clinic)

    #---------------------------------------------------------------------
    # results
    return df_Hosp, df_clinic
end
export HospClinic










############################ FUNCTIONS #############################
#----------------------------------------------
# reading infections and populations per agegroup
#-----------------------------------------------
function read_file(file::String)
    fname = file  #.csv file
    fpath = "/Users/shokoofehnourbakhsh/Dropbox/YorkUniversity/RSV_Shokoofeh/JuliaCodes/"
    path = fpath*fname
    f = CSV.read(path,normalizenames=true) #DataFrame table
    return f
end
export read_file


#-----------------------------------------------------------
# calculate incidence rate's denominator (person-days)/100,000
#-----------------------------------------------------------
function calcu_denominator(pop::Array)
    # population of 1)total,2)healthy,3)preterm_ill per agegroup i)0-2,ii)3-5,iii)6-11,iv)12-23,v)24-36 months

    days = 178 # period of one season in days; i.e. January to July
    denom = similar(pop,Float64) # data collection array
    denom = pop .* days/100000 #apply days/100000 on every pop
    #denom = mapcols(x-> x*days/100000 ,p) #apply days/100000 on every elements of p::DataFrame

    return denom
end
export calcu_denominator

#-------------------------------------------------------
# Calculate regional hosp. cases and outpatiant numbers
#-------------------------------------------------------
function calcu_HospClinic(S::Int64, inf::Array,denom::Array{Float64},inc::Array{Float64})
    # data-collection vectors
    Hosp = similar(inf) # Hosp cases ::Int64
    clinic = similar(inf) # clinic patients ::Int64

    #-----------------------------
    # calculate hospitalization cases (i.e. can not exceed maximum infected kids)
    _Hosp = similar(denom) # empty array for middle calculations (faster performance)
    _Hosp = denom .* inc # raw float number of hospital cases


    if S == 2
        # reduced risk intervals of hospitalizations by Palivizumab for preterm_ill children
        # agegroups (only,preterm_ill); 0-2mon. 3-5mon. 6-11mon. 12-23mon.)
        risk_lower = @SVector[0.05,0.05,0.05,0.05]
        risk_upper = @SVector[0.9626,0.9626,0.634,0.634]
        rows, cols = size(HospS1)
        for sim=1:rows
            Random.seed!(sim)
            for j in eachindex(risk_lower)
                risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                _Hosp[sim,j+j] = _h
                _h = _h * (1-risk)
            end
        end
    end

    Hosp = ceil.(Int64, _Hosp)
    @inbounds for i in eachindex(_Hosp)
        # ensure hospital admissions do not exceed the total infection
        Hosp[i] = minimum([_Hosp[i],inf[i]]) # hosp cases
    end
    #---------------------------
    # Whoever NOT admitted to the hospital is assumed to go to the local clinic
    # outpatient clinic cases
    clinic = inf - Hosp

    #---------------------------
    #output of the function
    return Hosp, clinic
end
export calcu_HospClinic



#-------------------------------------------------------------------------------
# Senario 2: INTERVENTION with Palivizumab for years before 2017 (old program).
# Prophylaxis are received only by 1- preterm infants < 6 months, 2- chronically ill kids < 2 years old
#-------------------------------------------------------------------

# inputs are regional hospital admission without any intervention (senario 1:S1)
function calcu_HospClinicS2(inf,HospS1)
    # data-collection vectors
    HospS2 = HospS1 # Hosp cases
    clinicS2 = similar(HospS1) # clinic patients


    # reduced risk intervals of hospitalizations by Palivizumab for preterm_ill children
    # agegroups (only,preterm_ill); 0-2mon. 3-5mon. 6-11mon. 12-23mon.)
    risk_lower = @SVector[0.05,0.05,0.05,0.05]
    risk_upper = @SVector[0.9626,0.9626,0.634,0.634]
    rows, cols = size(HospS1)
    for sim=1:rows
        Random.seed!(sim)
        for j in eachindex(risk_lower)
            risk = rand(Uniform(risk_lower[j],risk_upper[j]))
            _Hosp = HospS1[sim,j+j] * (1-risk)
            HospS2[sim,j+j] = ceil(_Hosp)
        end
    end

    #---------------------------
    # Whoever NOT admitted to the hospital is assumed to go to the local clinic
    # outpatient clinic cases
    clinicS2 = inf - HospS2

    #---------------------------
    #output of the function
    return HospS2, clinicS2

end
export calcu_HospClinicS2


#--------------------------------------
# to run in REPL
@everywhere using .RSV_Hosp
####################
end #module
