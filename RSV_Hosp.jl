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

#-----------------------------
# calculate hospitalization rate's denominator per simulation
# denom = person*days/100000
#denom = calcu_denominator(pop)


################# Main Function ##########################
# In order for this function to work we should have three files in the Julia working directory:
# infection.csv, incidence.csv, population.csv produced before by simulation's script
# Those files read by function read_file(); you need to customize fpath in read_file() to able to read files
function main_Hosp()
    denom = calcu_denominator()

    # calculate Regional Hosp, clinic's patients and
    # Tertiary Hosp that includs General Ward and ICU cases per all senarios
    # S1: No Paliv. S2: Paliv. > 2017, S3: Paliv. < 2017, S4: Maternal Vaccine
    # Note that LAMA produces same results as Palivizumab - no need for re-calculations
    for S=1:4
        Random.seed!(S)
        HospClinic_GWardICU(S,denom)
    end
end
export main_Hosp

############################ Used FUNCTIONS #############################

#-----------------------------------------------------------
# calculate incidence rate's denominator (person-days)/100,000
#-----------------------------------------------------------
function calcu_denominator()
    # population of 1)total,2)healthy,3)preterm_ill per agegroup i)0-2,ii)3-5,iii)6-11,iv)12-23,v)24-36 months
    # pop0to2_total | pop0to2_H | pop0to2_P | ... 500×13 DataFrame table
    population = read_file("population.csv") # population per simulation in  a row
    _pop = select(population, Not([1,4,7,10])) #remove columns (we only need healthy and preterm_ill population)
    pop = Matrix(_pop)

    # calculate hospitalization rate's denominator per simulation
    days = 178 # period of one season in days; i.e. January to July
    denom = similar(pop,Float64) # data collection array
    denom = (pop .* days)/100000 #apply days/100000 on every pop

    return denom
end
export calcu_denominator

#----------------------------------------------
# reading csv.files from saved-data of the simulation
#-----------------------------------------------
function read_file(file::String)
    fname = file  #.csv file
    fpath = "/Users/shokoofehnourbakhsh/Dropbox/YorkUniversity/RSV_Shokoofeh/JuliaCodes/"
    path = fpath*fname
    f = CSV.read(path,normalizenames=true) #DataFrame table
    return f
end
export read_file

#----------------------------------------------------------------------------
# Estimate Hospital Admissions and clinc outpatients for All senarios; 1,2,3,4
#-----------------------------------------------------------------------------
# S represent senarios;
# 1- No intervention 2- Palivizumab program < 2017, 3- Palivizumab program > 2017
# long-acting monoclonal antibody has the same efficacy as palivizumab (no need for seperate calculation)
# 4- Maternal vaccine (ResVax)
function HospClinic_GWardICU(S::Int64, denom::Array{Float64})

    #---------------------------
    # Reading total infection (results) from Agent-based Simulation
    # introduced_inf | inf0to2_total | inf0to2_H | inf0to2_P | ...500×14 DataFrame table
    infection = read_file("infection.csv") # infections per simulation in a row
    _inf = select(infection, Not([1,2,5,8,11])) # take out 'total infection' columns per agegroup
    inf = Matrix(_inf) # convert to matrix

    # introduced_inf | inc0to2_H | inc0to2_P | ...  500X10 DataFrame table
    incidence =  read_file("incidence.csv") # hosp. incidence rates per agegroup
    _inc = select(incidence, Not(1)) #take out initial infection column
    inc = Matrix(_inc)

    #------------------------------
    # Determine hospitalization cases and outpatient numbers per every simulation
    # S determines the senario
    Hospcases, outpatients = calcu_HospClinic(S,inf,denom,inc)

    #-------------------------
    # Tertiary Hospital admissions; general wards or ICU
    # S determines the senario
    ICU , GWard = calcu_GWardICU(S,Hospcases)


    #----------------------------------------------------------
    # DataFrame to collect data
    # names
    names_cases = Symbol.(["HospS$S 0to2_H", "HospS$S 0to2_P","HospS$S 3to5_H", "HospS$S 3to5_P",
    "HospS$S 6to11_H","HospS$S 6to11_P","HospS$S 12to23_H", "HospS$S 12to23_P", "HospS$S 24to35"])

    names_patients = Symbol.(["clinicS$S 0to2_H", "clinicS$S 0to2_P","clinicS$S 3to5_H", "clinicS$S 3to5_P",
    "clinicS$S 6to11_H","clinicS$S 6to11_P","clinicS$S 12to23_H", "clinicS$S 12to23_P", "clinicS$S 24to35"])

    names_GWard = Symbol.(["GWardS$S 0to2_H", "GWardS$S 0to2_P","GWardS$S 3to5_H", "GWardS$S 3to5_P",
    "GWardS$S 6to11_H","GWardS$S 6to11_P","GWardS$S 12to23_H", "GWardS$S 12to23_P", "GWardS$S 24to35"])

    names_ICU = Symbol.(["ICUS$S 0to2_H", "ICUS$S 0to2_P","ICUS$S 3to5_H", "ICUS$S 3to5_P",
    "ICUS$S 6to11_H","ICUS$S 6to11_P","ICUS$S 12to23_H", "ICUS$S 12to23_P", "ICUS$S 24to35"])

    df_Hosp = DataFrame(Hospcases, names_cases)
    df_clinic = DataFrame(outpatients, names_patients)
    df_GWard = DataFrame(GWard, names_GWard)
    df_ICU = DataFrame(ICU, names_ICU)

    #---------
    # saving data
    CSV.write("HospS$S.csv",df_Hosp)
    CSV.write("clinicS$S.csv",df_clinic)
    CSV.write("GWardS$S.csv",df_GWard)
    CSV.write("ICUS$S.csv",df_ICU)

    #---------------------------------------------------------------------
    # results
    return df_Hosp, df_clinic, df_GWard, df_ICU
end
export HospClinic_GWardICU


#-------------------------------------------------------------------------
# Calculate 1- Regional Hosp. and outpatiant cases plus
#---------------------------------------------------------------------------
function calcu_HospClinic(S::Int64, inf::Array{Int64},denom::Array{Float64},inc::Array{Float64})
    # data-collection vectors
    Hosp = similar(inf)   #Regional Hosp cases ::Int64
    clinic = similar(inf) #clinic patients ::Int64
    #-----------------------------
    # calculate hospitalization cases (i.e. can not exceed maximum infected kids)
    _Hosp = similar(denom) # empty array for middle calculations (faster performance)
    _Hosp = denom .* inc # raw Float numbers of hospital cases
    nrows, ncols = size(_Hosp)

    #-------------------------------------------------------------------------------
    # Senario 2,3: INTERVENTION with Palivizumab for years before(S2)/after(S3) 2017.
    # INTERVENTION with LAMA for years BEFORE(S4)/AFTER(S5) 2017 is the same as S2 and S3 (assuming same afficacy as Paliv.)
    # Before 2017: Prophylaxis are received only by 1- preterm infants < 6 months, 2- chronically ill kids < 2 years old
    # After 2017: Prophylaxis applied for both eligible infants 1 and 2 plus 3-healthy infants < 3 months
    #-------------------------------------------------------------------
    if S == 2 || S == 3
        # reduced risk intervals of hospitalizations by Palivizumab for preterm_ill children
        # agegroups (only preterm_ill); 0-2mon. 3-5mon. 6-11mon. 12-23mon.)
        risk_lower = @SVector[0.05,0.05,0.05,0.05]
        risk_upper = @SVector[0.9626,0.9626,0.634,0.634]
        for sim=1:nrows
            for j in eachindex(risk_lower)
                risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                _Hosp[sim,j+j] = _Hosp[sim,j+j] * (1-risk)
            end
        end

        #-------------------------------------------------------------------------------
        # Senario 3: INTERVENTION with Palivizumab for years AFTER 2017 (updated program).
        # Prophylaxis are received by:
        # 1- preterm infants < 6 months, 2- chronically ill kids < 2 years old
        # 3- healthy term infants < 3 months old
        #-------------------------------------------------------------------
        if S == 3
            for sim=1:nrows
                risk = rand(Uniform(0.05,0.50)) #reduced-risk factor interval for healthy terms
                _Hosp[sim,1] = _Hosp[sim,1] * (1-risk) #column 1 is healthy term 0-2months old
            end
        end


    #-------------------------------------------------------------------------------
    # Senario 4: INTERVENTION with ResVax (Maternal Vaccine).
    # Prophylaxis are received only by agegroup 0-2 months (both healthy and preterm_ill)
    #-------------------------------------------------------------------
    elseif S == 4
        # reduced risk intervals of hospitalizations through 90 days serveillance
        # agegroup 0-2mon. both healthy and preterm_ill
        risk_lower = @SVector[0.196,0.196]
        risk_upper = @SVector[0.615,0.615]
        for sim=1:nrows
            for j in eachindex(risk_lower)
                risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                _Hosp[sim,j] = _Hosp[sim,j] * (1-risk)
            end
        end
    end # senarios if

    #-------------------------
    # round it up to an integer value
    Hosp = ceil.(Int64, _Hosp)
    @inbounds for i in eachindex(Hosp)
        # ensure hospital admissions do not exceed the total infection
        Hosp[i] = minimum([Hosp[i],inf[i]]) # hosp cases
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


#-------------------------------------------------------------------
# calculate General Ward and ICU admissions for different senarios
#-------------------------------------------------------------------
function calcu_GWardICU(S::Int64,Hosp::Array)
    # collecting arrays
    ICU = similar(Hosp)
    TerHosp = similar(Hosp)
    GWard = similar(Hosp)
    _TerHosp = similar(Hosp,Float64)
    _ICU = similar(Hosp,Float64)

    nrows, ncols = size(Hosp)
    if S == 1 || S == 4
        # parameters
        TerHosp_rate = 0.955
        ICU_rate = 0.50

        # Tertiary hospital admissions
        _TerHosp = Hosp*TerHosp_rate
        #------------------------------------------
        # Maternal Vaccine (ResVax)
        if S == 4
            for sim=1:nrows
                risk = rand(Uniform(0.321,0.760)) #reduced-risk factor for both healthy/preterm_ill
                _TerHosp[sim,1] = _TerHosp[sim,1] * (1-risk) #column 1 is healthy term 0-2months old
                _TerHosp[sim,2] = _TerHosp[sim,2] * (1-risk) #column 2 is preterm_ill 0-2months old
            end
        end

    elseif S == 2 || S == 3
        # parameters
        ICU_rate = 0.50
        TerHosp_low = @SVector[0.22222,0.22222,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        TerHosp_up = @SVector[0.5,0.5,0.5,0.5,0.25,0.25,0.0001,0.0001,0.0001]
        rows, cols = size(Hosp)
        for sim=1:nrows
            for j in eachindex(TerHosp_low)
                rate = rand(Uniform(TerHosp_low[j],TerHosp_up[j]))
                _TerHosp[sim,j] = Hosp[sim,j] * rate
            end
        end

        #------------------------------------------
        # Paliv. > 2017; relative reduced risk by palivizumab for healthy terms in 0-2 months age
        if S == 3
            for sim=1:nrows
                risk = rand(Uniform(0.476,0.767)) #reduced-risk factor interval for healthy terms
                _TerHosp[sim,1] = _TerHosp[sim,1] * (1-risk) #column 1 is healthy term 0-2months old
            end
        end

    end #if

    # ICU Cases
    _ICU = _TerHosp*ICU_rate

    # round up ICU and Tertiary Hosp to integer
    ICU = ceil.(Int64,_ICU)
    TerHosp = ceil.(Int64,_TerHosp)

    # General ward admissions
    GWard = TerHosp - ICU

    #---------------------------
    #output of the function
    return ICU, GWard
end
export calcu_GWardICU



#--------------------------------------
# to run in REPL
@everywhere using .RSV_Hosp
####################
end #module
