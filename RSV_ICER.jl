## Cost-effectiveness for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh

module RSV_ICER

## Loading necessary packages

#using ClusterManagers
using Distributed
using Distributions #sample() lives here
#using StatsBase
using StaticArrays
#using Parameters
using DataFrames
using CSV
using Statistics
using Random
using Plots





################# Main Function #######################
# Main function to calculate ICER per simulation per agegroup
# for different baseline and intervention senarios

function main_ICER()
    senario_0 = @SVector[1,2,3]
    senario_1 = @SVector[1,2,3,4,5,6,7]
    l_0 = length(senario_0)
    l_1 = length(senario_1)

    # table to collect ICERs with different baseline senario of S1/S2/S3
    table = zeros(Float64,(l_0,l_1))
    CI_low = zeros(Float64,(l_0,l_1))
    CI_high = zeros(Float64,(l_0,l_1))
    @inbounds for S_0 = 1:l_0
        for S_1 = S_0+1:l_1
            ICER_agegroup, ICER_total = get_ICER(senario_0[S_0],senario_1[S_1]) #ICER with different bases

            boot_mean, CI_l, CI_h = bootstrap_mean(ICER_total)

            # average of total ICER per hospital stay using different intervention
            table[S_0,S_1] = mean(boot_mean)
            CI_low[S_0,S_1] = CI_l
            CI_high[S_0,S_1] = CI_h

            #----------------------------------------------------------
            # DataFrame to collect data
            # names
            names_age = Symbol.(["ICER_S$S_0$S_1,0to2_H", "ICER_S$S_0$S_1 0to2_P","ICER_S$S_0$S_1 3to5_H", "ICER_S$S_0$S_1 3to5_P",
            "ICER_S$S_0$S_1 6to11_H","ICER_S$S_0$S_1,6to11_P","ICER_S$S_0$S_1 12to23_H", "ICER_S$S_0$S_1 12to23_P", "ICER_S$S_0$S_1 24to35"])

            df_ICER_agegroup = DataFrame(ICER_agegroup,names_age)
            #df_ICER_total = DataFrame(ICER = ICER_total)

            #save ICER data
            CSV.write("ICER_agegroup$S_0$S_1.csv",df_ICER_agegroup)
            #CSV.write("ICER_total$S_0$S_1.csv",df_ICER_total)
        end
    end

    # save ICER table
    names_table = Symbol.(["S1:NoInter","S2:Paliv<2017","S3:Palivi>2017","S4:ResVax","S5:ResVax/Palivi","S6:LAMA<2017","S7LAMA>2017"])
    df_ICER_table = DataFrame(table, names_table)
    insertcols!(df_ICER_table, 1, baseline = ["S1","S2","S3"])
    CSV.write("ICER_table.csv",df_ICER_table)


    return df_ICER_table, CI_low, CI_high
end
export main_ICER


################# Used Functions #######################
#------------------------------------------------------------------------------
# Calculate ICER with respect to the S_0:baseline senario and S_1: campared senario
# ICER = (costS_1 - costS_0)/(effS_1 - effS_0)
#------------------------------------------------------------------------------
function get_ICER(S_0::Int64,S_1::Int64)

    #-------------- ICER for different agegroups and different simulation
    C1, E1 = CostEff(S_1)
    C0, E0 = CostEff(S_0)
    # calculate cost-effectiveness per hospital stay
    CE  = similar(C1) # data collecting array::Float64
    deltaC = C0 - C1
    deltaE = -(E0 - E1)
    CE  = deltaC ./ deltaE # for different agegroups per simulation
    #### (note that CE per agegroup contains NaN and Inf elements- I'll leave it to clean later)

    nrows, ncolm = size(CE)
    #-------------- Total ICER for different simulations
    C1_total , E1_total = CostEff_total(S_1)
    C0_total , E0_total = CostEff_total(S_0)
    # calculate ICER per simulation level
    CE_total = zeros(Float64,nrows)
    CE_total = (C0_total - C1_total) ./ -(E0_total - E1_total)

    #clean data from NaN and Inf (zeros in costs or hospital stays)
    deleteat!(CE_total, CE_total.== Inf)
    deleteat!(CE_total, CE_total.== -Inf)
    deleteat!(CE_total, CE_total.== NaN)

    return CE, CE_total
end
export get_ICER


#---------------------------------------------------------------------------
# Calculate total costs and total hospital days for different agegroups
#---------------------------------------------------------------------------
function CostEff(S::Int64)
    # calculate total costs per agegroup
    clinic = cost_clinic(S)
    Hosp, ICU, Hosp_days, ICU_days = cost_HospICU(S)

    if S == 1
        cost_total = clinic + Hosp + ICU
    else
        prophylaxis = cost_palivi_LAMA_ResVax(S)
        cost_total = clinic + Hosp +  ICU + prophylaxis
    end

    # calculate total hospital lengths of stay
    days_total = Hosp_days + ICU_days

    # results
    return cost_total, days_total
end
export CostEff



#---------------------------------------------------
# Calculate total costs/days per simulation level
#---------------------------------------------------
function CostEff_total(S::Int64)
    # get costs and effects for all agegroups
    C_age , E_age = CostEff(S)
    nrows, ncolm = size(C_age)
    C_total = [sum(C_age[sim,:]) for sim=1:nrows]
    E_total = [sum(E_age[sim,:]) for sim=1:nrows]

    return C_total, E_total
end
export CostEff_total



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

#----------------------------------------------------
# Calculate costs associate with patients in
# 1) clinic, 2)regional hosp + tertiary general ward, 3)tertiary ICU
#----------------------------------------------------
#------------ Local Clinic ----------------------
function cost_clinic(S::Int64)
    # LAMA and palivi have same efficacy and produces same Hosp/ICU results
    if S == 6 #LAMA < 2017
        S = 2
    elseif S == 7 #LAMA > 2017
        S = 3
    end

    # average cost per case for mild RSV infection (no hospital admission)
    cost_mild = 1569

    # reading data file from my working directory
    _clinic = read_file("clinicS$S.csv")
    clinic = Matrix(_clinic) # convert DataFrame to matrix

    # calculate total costs
    cost = clinic .* cost_mild
    return cost
end
export cost_clinic


#------------------ Regional Hospital ---------------
function cost_HospICU(S::Int64)
    # LAMA and palivi have same efficacy and produces same Hosp/ICU results
    if S == 6 # LAMA < 2017
        S = 2
    elseif S == 7 #LAMA > 2017
        S = 3
    end
    #------- Parameters
    # ------- AVERAGE COSTS per inpatient case ended up to hospital or ICU
    cost_Hosp = 16946
    cost_ICU = 80590

    # range for minimum stay in hospital or ICU per agegroups
    HospStay = @SVector[5:8,9:19,5:8,9:19,5:8,9:19,  0:0,0:0,0:0]
    ICUstay = @SVector[5:9,11:16,5:9,11:16,5:9,11:16,   0:0,0:0,0:0]

    #-------- reading data file from my working directory
    _Hosp = read_file("HospS$S.csv")
    Hosp = Matrix(_Hosp) # convert DataFrame to Matrix

    _ICU = read_file("ICUS$S.csv")
    ICU = Matrix(_ICU) # convert DataFrame to Matrix



    #--------------------- total length of hospital days per agegroup
    # hospital stays per individual are selected and added up for each agegroup
    daysHosp = similar(Hosp,Int64) #collecting array
    nrows, ncolm = size(Hosp)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            daysHosp[sim,age] = sum(rand(HospStay[age],Hosp[sim,age]))
        end
    end
    # ICU stays per individual are selected and added up for each agegroup
    daysICU = similar(ICU,Int64) # collecting array
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            daysICU[sim,age] = sum(rand(ICUstay[age],ICU[sim,age]))
        end
    end


    #-------- total cost for all hospital/ICU admissions per agegroup
    # note: Hosp includes patients that might end up in ICU - their costs consists in cost_ICU
    costsHosp = similar(Hosp,Float64) #collecting array
    costsHosp = (Hosp - ICU) .* cost_Hosp # medium-infected cases

    costsICU = similar(ICU,Float64) #collecting array
    costsICU = ICU .* cost_ICU # sever infected cases

    return costsHosp, costsICU, daysHosp, daysICU
end
export cost_HospICU




#-----------------------------------------------------------
# Costs associate with intervention of palivizumab
#-----------------------------------------------------------
function cost_palivi_LAMA_ResVax(S::Int64)
    #------ parameters
    # Palivizumab costs per dose for different agegroup and conditions
    # average cost per dose plus $50 administration fee
    paliviCost_p = @SVector[1065.32,1567,2048.4,0.0] # preterm_ill for agegroups: 0-2,3-5,6-11,12-23
    #later should replace zero for agegroup 12-23=2657.4
    paliviCost_H = 1065.33 # healthy term for agegroup 0-2 months.

    # assumption cost for maternal vaccine (ResVax)
    ResVax_p = @SVector[1065.32,1567,2048.4,0.0]
    ResVax_H = 1065.33

    #------ read data file population from my woring directory
    population = read_file("population.csv")
    _popu = select(population, Not([1,4,7,10])) # take out 'total population' columns per agegroup
    popu = Matrix(_popu)

    l_p = length(paliviCost_p)
    total = zeros(Float64,size(popu)) # collecting array
    nrow, ncolm = size(total)

    # Maternal Vaccine (ResVax)
    if S ==4 || S == 5
        S == 4 #ResVax only
        @inbounds for sim=1:nrow
            Random.seed!(sim)
            total[sim,1] = ResVax_H # healthy term < 3 months
            total[sim,2] = ResVax_p[1] # preterm < 3 months
        end

        # ResVax + Palivi
        if S == 5
            # five doses for only preterm_ill > 3 months
            @inbounds for sim=1:nrow
                Random.seed!(sim)
                for j=2:l_p
                    total[sim,j+j] = 5 * paliviCost_p[j]
                end
            end
        end

    # palivi. or LAMA < 2017
    elseif S == 2 || S == 6
        @inbounds for j=1:l_p
            for sim=1:nrow
                Random.seed!(sim)
                # five doses for only preterm_ill
                # no cost for healthy terms; remained zero
                if S == 2 #Palivi
                    total[sim,j+j] = 5 .* paliviCost_p[j]
                elseif S == 6 #LAMA
                    total[sim,j+j] = paliviCost_p[j]
                end
            end
        end

    # palivi. or LAMA > 2017
    elseif S == 3 || S == 7
        @inbounds for sim=1:nrow
            Random.seed!(sim)
            if S == 3 #Palivi
                # cost for healthy terms < 3 months old
                total[sim,1] = 5 .* paliviCost_H
                for j=1:l_p
                    # five doses for only preterm_ill
                    total[sim,j+j] = 5 .* paliviCost_p[j]
                end
            elseif S == 7 #LAMA
                # cost for healthy terms < 3 months old
                total[sim,1] = paliviCost_H
                for j=1:l_p
                    # five doses for only preterm_ill
                    total[sim,j+j] = paliviCost_p[j]
                end
            end
        end
    end

    # total prophylaxis costs per simulation for different agegroups
    costs = similar(total) # collecting array
    costs = popu .* total

    return costs
end
export cost_palivi_LAMA_ResVax



#---------------------------------------------------------------------
# Bootstrap my ICER results from 500 simulation to get mean distribution
#---------------------------------------------------------------------
function bootstrap_mean(ICER::Array{Float64,1})
    leng = length(ICER)
    ICER_mean = mean(ICER)
    n_boot = 200
    boot_mean = zeros(Float64,n_boot)
    boot_delta = similar(boot_mean)
    @inbounds for i=1:n_boot
        # mean of ith bootstrap sample
        m = mean(sample(ICER,leng; replace=true))
        boot_mean[i] = m

        # CI95% - confidence interval
        boot_delta[i] = m - ICER_mean # boot_mean minus average of actual sample
    end

    # CI:95% - confidence interval
    s = sort(boot_delta)
    delta_head = s[Int(n_boot * 0.975)]
    delta_tail = s[Int(n_boot * 0.025)]

    CI_low = ICER_mean - delta_head
    CI_high = ICER_mean - delta_tail

    return boot_mean, CI_low, CI_high
end
export bootstrap_mean




#--------------------------------------------------------
# Histogram Plot of BOOTSTRAP mean distribution of ICERS
#-------------------------------------------------------
function histo(ICER)
    boot_mean, CI_low, CI_high = bootstrap_mean(ICER)
    boot_ave = mean(boot_mean)
    h = histogram(boot_mean,
                    title = "Bootstrap Distribution of median ICERs",
                    xlabel = "mean",
                    ylabel = "Frequency",
                    legend = false
                    )
    vline!(h, [CI_low,CI_high], line = (:green,2))
    vline!(h,[boot_ave], line = (:black,2))
end
export histo





#--------------------------------------
# to run in REPL
@everywhere using .RSV_ICER
####################
end #module
