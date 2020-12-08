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
#using Plots
using Bootstrap
using StatsBase





################# Main Function #######################
# Main function to calculate ICER per simulation per agegroup
# for different baseline and intervention senarios

function main_ICER()
    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        senario_0 = @SVector[1,2,3]
        senario_1 = @SVector[1,2,3,4,5,6,7]
        l_0 = length(senario_0)
        l_1 = length(senario_1)

        # table to collect ICERs with different baseline senario of S1/S2/S3
        table = zeros(Float64,(l_0,l_1))
        CI_min = similar(table,Float64)
        CI_max = similar(table,Float64)
        @inbounds for S_0 in senario_0
            for S_1 = S_0+1:l_1
                if S_0==2 && S_1==6 continue end
                if S_0==3 && S_1==7 continue end

                ICER_total = get_ICER(S_0,S_1,season) #ICER with different bases

                # boot_mean, CI_l, CI_h = bootstrap_mean(ICER_total)
                # average of total ICER per hospital stay per intervention
                b = bootstrap(mean,ICER_total,BasicSampling(2000))
                CI = confint(b, BCaConfInt(0.95)) #bias-corrected and accelerated CI
                table[S_0,S_1] = CI[1][1] #mean
                CI_min[S_0,S_1] = CI[1][2] # 95% CI_min
                CI_max[S_0,S_1] = CI[1][3] # 95% CI_max
            end #base scenario
        end #alternative scenario

        # save ICER table
        names_table = Symbol.(["S0_mean","S1_mean","S2_mean","S3_mean","S4_mean","S5_mean","S6_mean"])
        names_min = Symbol.(["S0_min","S1_min","S2_min","S3_min","S4_min","S5_min","S6_min"])
        names_max = Symbol.(["S0_max","S1_max","S2_max","S3_max","S4_max","S5_max","S6_max"])
        df_ICER_table = DataFrame(table,names_table)
        df_ICER_CImin = DataFrame(CI_min,names_min)
        df_ICER_CImax = DataFrame(CI_max,names_max)
        insertcols!(df_ICER_table, 1, :baseline => ["S0","S1","S2"], makeunique=false)
        insertcols!(df_ICER_CImin, 1, :baseline => ["S0","S1","S2"], makeunique=false)
        insertcols!(df_ICER_CImax, 1, :baseline => ["S0","S1","S2"], makeunique=false)
        CSV.write("ICER_table_$season.csv",df_ICER_table)
        CSV.write("ICER_CImin_$season.csv",df_ICER_CImin)
        CSV.write("ICER_CImax_$season.csv",df_ICER_CImax)
    end #season
end
export main_ICER



################# Used Functions #######################
#------------------------------------------------------------------------------
# Calculate ICER with respect to the S_0:baseline senario and S_1: campared senario
# ICER = -(costS0 - costS1)/(effS0 - effS1)
#------------------------------------------------------------------------------
function get_ICER(S_0::Int64,S_1::Int64,season::Symbol)
    #-------------- Total ICER per simulation
    C1_total, E1_total = CostEff_total(S_1,season) #costs and effects for S1
    C0_total, E0_total = CostEff_total(S_0,season) #costs and effects for base (S0)
    # calculate ICER per simulation level
    CE_total = (C0_total - C1_total) ./ -(E0_total - E1_total)

    #clean data from NaN and Inf (zeros in costs or hospital stays)
    deleteat!(CE_total, CE_total.== Inf)
    deleteat!(CE_total, CE_total.== -Inf)
    deleteat!(CE_total, CE_total.== NaN)

    #save raw ICER values per simulation for plotting
    df_ICER = DataFrame(ICERs=CE_total)
    CSV.write("ICER_S$S_0$S_1$season.csv",df_ICER)
    return CE_total
end
export get_ICER


function incremental()
    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        senario_0 = @SVector[1,2,3]
        senario_1 = @SVector[1,2,3,4,5,6,7]
        l_0 = length(senario_0)
        l_1 = length(senario_1)

        #collecting file
        inc_cost = zeros(500,15)
        inc_day = similar(inc_cost)
        i = 1
        @inbounds for S_0 in senario_0
            for S_1 = S_0+1:l_1
                inc_cost[:,i], inc_day[:,i] = get_incremental(S_0,S_1,season)
                i += 1
            end # baseline
        end #alternative scenario

        names = Symbol.(["S0S1","S0S2","S0S3","S0S4","S0S5","S0S6",
        "S1S2","S1S3","S1S4","S1S5","S1S6",
        "S2S3","S2S4","S2S5","S2S6"])
        df_cost = DataFrame(inc_cost,names)
        df_day = DataFrame(inc_day,names)
        CSV.write("inc_cost_$season.csv",df_cost)
        CSV.write("inc_day_$season.csv",df_day)

    end #season
end
export incremental


function get_incremental(S_0::Int64,S_1::Int64,season::Symbol)
    #-------------- Total incremental per simulation
    C1_total, E1_total = CostEff_total(S_1,season) #costs and effects for S1
    C0_total, E0_total = CostEff_total(S_0,season) #costs and effects for base (S0)
    # calculate incrementals per simulation level
    inc_costs = C1_total - C0_total
    inc_effects = E1_total - E0_total

    #save incremental values per simulation for cost-effectivenss plane
    #df_incremental = DataFrame(delta_costs=inc_cost, delta_days = inc_effect)
    #CSV.write("incremental_S$S_0$S_1$season.csv",df_incremental)
    return inc_costs, inc_effects
end
export get_incremental

#---------------------------------------------------------------------------
# Calculate total costs and total hospital days for different agegroups
#---------------------------------------------------------------------------
function CostEff(S::Int64,season::Symbol)
    # calculate total costs per agegroup
    clinic = cost_clinic(S,season)
    Hosp, ICU, Hosp_days, ICU_days = cost_HospICU(S,season)

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
function CostEff_total(S::Int64,season::Symbol)
    # get costs and effects for all agegroups and health status at the birth
    C_age , E_age = CostEff(S,season)
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
function cost_clinic(S::Int64,season::Symbol)
    # LAMA and palivi have same efficacy and produces same Hosp/ICU results
    if S == 6 #LAMA < 2017
        S = 2
    elseif S == 7 #LAMA > 2017
        S = 3
    end

    # average cost per case for mild RSV infection (no hospital admission)
    cost_mild = 1569

    # reading data file from my working directory
    _clinic = read_file("clinicS$S$season.csv")
    clinic = Matrix(_clinic) # convert DataFrame to matrix

    # calculate total costs
    cost = clinic .* cost_mild
    return cost
end
export cost_clinic


#------------------ Regional Hospital ---------------
function cost_HospICU(S::Int64, season::Symbol)
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
    _Hosp = read_file("HospS$S$season.csv")
    Hosp = Matrix(_Hosp) # convert DataFrame to Matrix

    _ICU = read_file("ICUS$S$season.csv")
    ICU = Matrix(_ICU) # convert DataFrame to Matrix



    #--------------------- total length of hospital days per agegroup
    # hospital and ICU stays per individual are selected and added up for each agegroup
    daysHosp = similar(Hosp,Int64) #collecting array
    daysICU = similar(ICU,Int64)
    nrows, ncolm = size(Hosp)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            daysHosp[sim,age] = sum(rand(HospStay[age],Hosp[sim,age]))
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
            total[sim,1] = 1560 # healthy term < 3 months
            total[sim,2] = 1560 # preterm < 3 months
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
    c = 0.95
    delta_head = s[Int(n_boot * (1+c)/2)]
    delta_tail = s[Int(n_boot * (1-c)/2)]

    CI_low = ICER_mean - delta_head
    CI_high = ICER_mean - delta_tail

    return boot_mean, CI_low, CI_high
end
export bootstrap_mean
"""
z = 1.959964
m = mean(ICER_mean)
s = stdm(ICER_mean, m; corrected=false)
m + z*(s/sqrt(length(ICER_mean)))
m - z*(s/sqrt(length(ICER_mean)))
"""
#--------------------------------------------------------
# Histogram Plot of BOOTSTRAP mean distribution of ICERS
#-------------------------------------------------------
function histo(ICER)
    boot_mean, CI_low, CI_high = bootstrap_mean(ICER)
    boot_ave = mean(boot_mean)
    h = histogram(boot_mean,
                    title = "Bootstrap Distribution of mean ICERs",
                    xlabel = "mean",
                    ylabel = "Frequency",
                    legend = false
                    )
    vline!(h, [CI_low,CI_high], line = (:green,2))
    vline!(h,[boot_ave], line = (:black,2))
end
export histo




#------------------------------------------------------------------------
# independent of the main ICER function save total costs and days for plotting
function get_costs_days()
    #collecting array
    c_total = zeros(Float64,(500,7))
    d_total = similar(c_total)

    sea = [:mild,:moderate,:severe]
    for season in sea
        for sc = 1:7
            costs,days = CostEff_total(sc, season)
            c_total[:,sc] = costs
            d_total[:,sc] = days
        end
        # saving the file
        names = Symbol.(["S$S" for S=1:7]) #[:S1,:S2..., :S7]
        df_c = DataFrame(c_total,names)
        df_d = DataFrame(d_total,names)
        CSV.write("costs$season.csv",df_c)
        CSV.write("days$season.csv",df_d)
    end
end
export get_costs_days
#------------------------------------------------
#--------------------------------------
# to run in REPL
@everywhere using .RSV_ICER
####################
end #module
