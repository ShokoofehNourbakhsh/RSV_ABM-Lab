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


################# Main Function ##########################
# In order for this function to work we should have three files in the Julia working directory:
# infection.csv, incidence.csv, population.csv produced before by simulation's script
# Those files read by function read_file(); you need to customize fpath in read_file() to able to read files
function main_Hosp()
    denom = calcu_denominator()

    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        # calculate Regional Hosp, clinic's patients and
        # Tertiary Hosp that includs General Ward and ICU cases per all senarios
        # S1: No Paliv. S2: Paliv. > 2017, S3: Paliv. < 2017, S4: Maternal Vaccine, S5: ResVax+Palivizumab
        # Note that LAMA produces same results as Palivizumab - no need for re-calculations
        for S=1:5
            ClinicHospICU(S,season,denom)
        end #scenarios
    end #season
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
    days = 180 # period of one season in days; i.e. January to July
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
function ClinicHospICU(S::Int64,season::Symbol,denom::Array{Float64})

    #---------------------------
    # Reading total infection (results) from Agent-based Simulation
    # introduced_inf | inf0to2_total | inf0to2_H | inf0to2_P | ...500×14 DataFrame table
    infection = read_file("infection_$season.csv") # infections per simulation in a row
    _inf = select(infection, Not([1,2,5,8,11])) # take out 'total infection' columns per agegroup
    inf = Matrix(_inf) # convert to matrix

    ######### due to lack of hospitalization data for 2 and 3 years old, put zeros on their respective columns
    nrows, ncolm = size(inf)
    inf[:,[9,8,7]] = zeros(Int64,(nrows,3))
    ################

    # introduced_inf | inc0to2_H | inc0to2_P | ...  500X10 DataFrame table
    incidence =  read_file("incidence_$season.csv") # hosp. incidence rates per agegroup
    _inc = select(incidence, Not(1)) #take out initial infection column
    inc = Matrix(_inc)

    #------------------------------
    # Determine hospitalization cases and outpatient numbers per every simulation
    # S determines the senario
    Hospcases, outpatients = calcu_HospClinic(S,inf,denom,inc)

    #-------------------------
    # Tertiary Hospital admissions which goe to ICU
    # S determines the senario
    ICU = calcu_ICU(S,Hospcases)


    #----------------------------------------------------------
    # DataFrame to collect data
    # names
    names_cases = Symbol.(["HospS$S 0to2_H", "HospS$S 0to2_P","HospS$S 3to5_H", "HospS$S 3to5_P",
    "HospS$S 6to11_H","HospS$S 6to11_P","HospS$S 12to23_H", "HospS$S 12to23_P", "HospS$S 24to35"])

    names_patients = Symbol.(["clinicS$S 0to2_H", "clinicS$S 0to2_P","clinicS$S 3to5_H", "clinicS$S 3to5_P",
    "clinicS$S 6to11_H","clinicS$S 6to11_P","clinicS$S 12to23_H", "clinicS$S 12to23_P", "clinicS$S 24to35"])

    names_ICU = Symbol.(["ICUS$S 0to2_H", "ICUS$S 0to2_P","ICUS$S 3to5_H", "ICUS$S 3to5_P",
    "ICUS$S 6to11_H","ICUS$S 6to11_P","ICUS$S 12to23_H", "ICUS$S 12to23_P", "ICUS$S 24to35"])

    df_Hosp = DataFrame(Hospcases, names_cases)
    df_clinic = DataFrame(outpatients, names_patients)
    df_ICU = DataFrame(ICU, names_ICU)

    #---------
    # saving data
    CSV.write("Hosp_S$S$season.csv",df_Hosp)
    CSV.write("clinic_S$S$season.csv",df_clinic)
    CSV.write("ICU_S$S$season.csv",df_ICU)

    #---------------------------------------------------------------------
    # results
    return df_clinic, df_Hosp, df_ICU
end
export ClinicHospICU


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
    _Hosp = denom .* inc # raw Float numbers of hospital cases (No Intervention)
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
        risk_lower = @SVector[0.20,0.20,0.20,0.20]
        risk_upper = @SVector[0.90,0.90,0.67,0.67]
        for sim=1:nrows # each simulation
            Random.seed!(sim)
            for j in eachindex(risk_lower) # each agegroup
                risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                _Hosp[sim,j+j] = _Hosp[sim,j+j] .* (1-risk)
            end
        end

        #-------------------------------------------------------------------------------
        # Senario 3: INTERVENTION with Palivizumab for years AFTER 2017 (updated program).
        # Prophylaxis are received by:
        # 1- preterm infants < 6 months, 2- chronically ill kids < 2 years old
        # 3- healthy full-term < 5 months old
        #-------------------------------------------------------------------
        if S == 3
            for sim=1:nrows
                # clinical reduced rate for healthy full-term between preintervention (2014-16) and intervention years (2017-19)
                # data obtained from https://www.sciencedirect.com/science/article/pii/S221133552030139X
                risk_1 = 0.235 # 0-2 months
                risk_3 = 0.38  # 3-5 months
                _Hosp[sim,1] = _Hosp[sim,1] .* (1-risk_1)
                _Hosp[sim,3] = _Hosp[sim,3] .* (1-risk_3)
            end
        end


    #-------------------------------------------------------------------------------
    # Senario 4: INTERVENTION with ResVax (Maternal Vaccine).
    # Prophylaxis are received only by agegroup 0-2 months (both healthy and preterm_ill)
    #-------------------------------------------------------------------
    elseif S == 4
        # reduced risk intervals of hospitalizations through 90 days serveillance
        # agegroup 0-2mon. both healthy and preterm_ill
        risk_lower = @SVector[0.247,0.319]
        risk_upper = @SVector[0.619,0.750]
        for sim=1:nrows
            Random.seed!(sim)
            for j in eachindex(risk_lower)
                risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                _Hosp[sim,j] = _Hosp[sim,j] .* (1-risk)
            end
        end

    #-------------------------------------------------------------------------------
    # Senario 5: INTERVENTION with ResVax + Palivi.
    # 1- ResVax protects 0-2 months (both healthy and preterm_ill)
    # 2- palivizumab shots will be provided for preterm_ill > 3 months but below 1 year-old
    #-------------------------------------------------------------------
    elseif S == 5
            # reduced risk intervals of hospitalizations through 90 days serveillance
            # agegroup 0-2mon. both healthy and preterm_ill
            risk_lower = @SVector[0.247,0.319]
            risk_upper = @SVector[0.619,0.750]
            for sim=1:nrows
                Random.seed!(sim)
                for j in eachindex(risk_lower)
                    risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                    _Hosp[sim,j] = _Hosp[sim,j] .* (1-risk)
                end
            end

            # reduced risk intervals of hospitalizations by Palivizumab for preterm_ill children
            # agegroups (only preterm_ill); 3-5mon. 6-11mon. 12-23mon.)
            risk_lower = @SVector[0.20,0.20,0.20]
            risk_upper = @SVector[0.90,0.67,0.67]
            for sim=1:nrows # each simulation
                Random.seed!(sim)
                for j in eachindex(risk_lower) # each agegroup
                    risk = rand(Uniform(risk_lower[j],risk_upper[j]))
                    _Hosp[sim,j+j+2] = _Hosp[sim,j+j+2] * (1-risk)
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
# calculate ICU admissions for different senarios
#-------------------------------------------------------------------
function calcu_ICU(S::Int64,Hosp::Array)
    # collecting arrays
    ICU = similar(Hosp)
    _ICU = similar(Hosp,Float64)
    _ICUResVax = similar(Hosp,Float64)

    # S1: No Intervention
    # portion of ICU admission out of hospitalized cases
    ICU_low = @SVector[0.06,0.0,  0.0,0.0,           0.0,0.0,   0.0,0.0,0.0]
    ICU_up = @SVector[0.25,1.0,   0.00001,0.00001,   0.17,1.0,   0.000001,0.0000001,0.0000001]
    nrows, ncols = size(Hosp)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for j in eachindex(ICU_low)
            rate = rand(Uniform(ICU_low[j],ICU_up[j]))
            _ICU[sim,j] = Hosp[sim,j] * rate
        end
    end

    #------------------------------------
    # Palivizumab before/after 2017
    if S == 2 || S == 3
        risk_p = 0.6390 # reduced risk preterm_ill < 12 months old
        risk_H = 0.4385 # reduced risk HFT < 12 months old
        @inbounds for j=1:3 # only preterm_ill < 12 months
            _ICU[:,[j+j]] = _ICU[:,j+j] .* (1-risk_p)
        end
        #-------------
        # Paliv. > 2017;
        if S == 3
            _ICU[:,1] = _ICU[:,1] .* (1-risk_H) # healthy term 0-2 months old
            _ICU[:,3] = _ICU[:,3] .* (1-risk_H) # healthy term 3-5 months old
        end

    #-------------------------------------------
    # Maternal Vaccine (ResVax) or ResVax+Palivi.
    elseif S == 4 || S == 5
        @inbounds for sim=1:nrows
            Random.seed!(sim)
            risk_H, risk_p = rand(Uniform(0.319,0.750),2) #reduced-risk factor for both healthy/preterm_ill
            _ICU[sim,1] = _ICU[sim,1] * (1-risk_H) # healthy term 0-2months old
            _ICU[sim,2] = _ICU[sim,2] * (1-risk_p) # preterm_ill 0-2months old
        end

        #------------
        # ResVax + Palivi
        if S == 5
            risk_p = 0.6390 # reduced risk preterm_ill < 12 months old
            @inbounds for j=2:3 # only preterm_ill between 3-5 and 6-11 months old
                _ICU[:,[j+j]] = _ICU[:,j+j] .* (1-risk_p)
            end
        end
    end #if

    #------------------------
    # round up ICU to integer
    ICU = round.(Int64,_ICU)

    #---------------------------
    #output of the function
    return ICU
end
export calcu_ICU



#--------------------------------------
# to run in REPL
@everywhere using .RSV_Hosp
####################
end #module
