## Cost-effectiveness for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh

module RSV_CEA

## Loading necessary packages

using Distributed
using Distributions #sample()
#using StatsBase
using StaticArrays
using DataFrames
using CSV
using Statistics
using Random
#using Plots
using Bootstrap
using StatsBase


################# Main Function #######################
#------------------------------------------------------------------------------
# Calculate ICER with respect to the S_0:baseline senario and S_1: campared senario
# ICER = -deltaCost/deltaEffect
#------------------------------------------------------------------------------

function bootdata_increments_ICER()
    # We will bootstrap simulation results of costs and days, then obtain their difference and then difference ratio (ICER)
    #bootstrapping data by season
    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        #load simulation data as a table
        _dat = read_file("CostEff_datatable_$season.csv")
        _dat = select(_dat,Not(15:35)) # keep only total costs and total QALY
        dat = Matrix(_dat)

        # prepare data for bootsrap function to iterate by coloumn
        # we pass a table of 500 original simulations of costs and effects
        # and perform bootstrap simultaniously for all costs/effects pairs
        r,c = size(dat)
        d = [dat[i,:] for i=1:r] # put all columns of one simulation in one vector

        # mean bootstrapping costs and days pairs for all scenarios
        Random.seed!(123)
        n_boot = 1000
        b = bootstrap(bootfc_increments_ICER,d,BasicSampling(n_boot))

        # bootfc_increments_ICER produces vector of 3*6 (6 is total compared scenarios)
        # seperate results from b into ICER, delta_cost, delta_day
        boot_ICER = zeros(Float64,(n_boot,6))
        boot_deltaCost = similar(boot_ICER)
        boot_deltaQ = similar(boot_ICER)
        # table for reporting final ICER result plus 95% confidence interval
        ICER_table = zeros(Int64,(3,6))
        for j=1:6
            boot_ICER[:,j] = b.t1[j]
            boot_deltaCost[:,j] = b.t1[j+6]
            boot_deltaQ[:,j] = b.t1[j+12]

            # we only want 95% CI of ICER (choosing index 1:15)
            # outputs of confint() are mean, lower_CI and uppper bound CI
            # BCaConfInt() means bias-corrected accelerated method
            mmean, CIlow, CIup = confint(b, BCaConfInt(0.95),j)
            ICER_table[1,j] = ceil(Int64,mmean)
            ICER_table[2,j] = ceil(Int64,CIlow)
            ICER_table[3,j] = ceil(Int64,CIup)

        end

        # save bootstrapped data
        names = Symbol.(["S0S1","S0S2","S0S3","S1S4","S2S5","S3S6"])
        df_bootICER = DataFrame(boot_ICER,names)
        df_bootDeltaCost = DataFrame(boot_deltaCost,names)
        df_bootDeltaQ = DataFrame(boot_deltaQ,names)
        df_ICERtable = DataFrame(ICER_table,names)

        CSV.write("bootICER_$season.csv",df_bootICER)
        CSV.write("bootDeltaCost_$season.csv",df_bootDeltaCost)
        CSV.write("bootDeltaQALY_$season.csv",df_bootDeltaQ)
        CSV.write("ICERtable_$season.csv",df_ICERtable)
    end
end
export bootdata_increments_ICER

# ------ function used for bootstrap() in above
# statistic of the interest
# bootfc will perform mean, difference of means (alternative - base),ratio of differences (delta_cost/delta_QALY)
function bootfc_increments_ICER(d)
    m = mean(d)
    scenario_0 = @SVector[1]
    scenario_1 = @SVector[2,3,4]

    #collecting file
    cost_inc = zeros(Float64,6)
    Q_inc = similar(cost_inc)
    ICER = similar(cost_inc)
    i = 1
    @inbounds for S_0 in scenario_0
        for S_1 in scenario_1
            delta_c = m[S_1] - m[S_0]
            delta_q = m[S_1+7] - m[S_0+7]
            cost_inc[i] = delta_c
            Q_inc[i] = delta_q
            ICER[i] = delta_c / delta_q
            i += 1
        end # baseline
    end #alternative

    S_0 = 2
    S_1 = 5
    delta_c = m[S_1] - m[S_0]
    delta_q = m[S_1+7] - m[S_0+7]
    cost_inc[i] = delta_c
    Q_inc[i] = delta_q
    ICER[i] = delta_c / delta_q
    i += 1

    S_0 = 3
    S_1 = 6
    delta_c = m[S_1] - m[S_0]
    delta_q = m[S_1+7] - m[S_0+7]
    cost_inc[i] = delta_c
    Q_inc[i] = delta_q
    ICER[i] = delta_c / delta_q
    i += 1

    S_0 = 4
    S_1 = 7
    delta_c = m[S_1] - m[S_0]
    delta_q = m[S_1+7] - m[S_0+7]
    cost_inc[i] = delta_c
    Q_inc[i] = delta_q
    ICER[i] = delta_c / delta_q

    # at the end, bootfc delivers big vector of 6 points of ICERs, 6 points of delta_cost and 6 points of delta_q
    return append!(ICER, append!(cost_inc,Q_inc))
end
export bootfc_increments_ICER




"""
###### COST MINIMIZED -------
function bootdata_costminimized()
    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        #load simulation data as a table
        _dat = read_file("CostEff_S1256season.csv")
        dat = Matrix(_dat)
        r,c = size(dat)
        d = [dat[i,:] for i=1:r] # put all columns of one simulation in one vector

        # mean bootstrapping costs for cost_minimized scenarios
        Random.seed!(123)
        n_boot = 1000
        b = bootstrap(bootfc_costminimized,d,BasicSampling(n_boot))

        # table for reporting final ICER result plus 95% confidence interval
        CM_table = zeros(Int64,(3,2))
        for j=1:2
            mmean, CIlow, CIup = confint(b, BCaConfInt(0.95),j)
            CM_table[1,j] = ceil(Int64,mmean)
            CM_table[2,j] = ceil(Int64,CIlow)
            CM_table[3,j] = ceil(Int64,CIup)
        end

        # save bootstrapped data
        names = Symbol.(["S1S5","S2S6"])
        df_CMtable = DataFrame(CM_table,names)
        CSV.write("CMtable_season.csv",df_CMtable)
    end
end
export bootdata_costminimized

# boot function for delta_cost in cost-minimized
function bootfc_costminimized(d)
    m = mean(d)

    cost_inc = zeros(Float64,2)
    # function
    cost_inc[1] = m[3] - m[1] #S1S5
    cost_inc[2] = m[4] - m[2] #S2S6

    # at the end, bootfc delivers big vector of 2 points of delta_costS1S5 and delta_costS2S6
    return cost_inc
end
export bootfc_increments_ICER

function bfc(d)
    m=mean(d)
    delta_cost = m[2]-m[1]
    return delta_cost
end
export bfc

function bootbfc()
    #_dat = read_file("CostEff_S1S3mild.csv")
    _dat = read_file("CostEff_S2S4moderate.csv")

    dat = Matrix(_dat)
    r,c = size(dat)
    d = [dat[i,:] for i=1:r] # put all columns of one simulation in one vector

    # mean bootstrapping costs for cost_minimized scenarios
    Random.seed!(123)
    n_boot = 1000
    b = bootstrap(bfc,d,BasicSampling(n_boot))
    mmean, CIlow, CIup = confint(b, BCaConfInt(0.95))

    return mmean,CIlow, CIup

end
export bootbfc
"""
################# Used Functions ########################------------------------------------------------------------------
#-------------- create data table of raw costs and days per scenarios (simulation raw data)
# MAIN function to create table of total costs and total QALY per scenario
function CostEff_datatable()
    sea = [:mild,:moderate,:severe]
    for season in sea

        inf = get_inf(season)

        #collecting array (costs for 7 scenarios and days for 7 scenarios)
        #table = zeros(Float64,(500,14))
        table = zeros(Float64,(500,35))
        for S= 1:7
            costs, QALY, Q_clinic, Q_GW, Q_ICU  = CostEff_total(S, season)
            if S == 1
                Q_averted = 0
            else
                Q_averted = Eff_averted(S,season,inf)
            end
            table[:,S]   = costs
            table[:,S+7] = QALY .+ Q_averted
            table[:,S+14]= Q_clinic
            table[:,S+21]= Q_GW
            table[:,S+28]= Q_ICU
        end
        # saving the file
        names = Symbol.(["costS$S" for S=1:7]) #[:costS1,:costS2..., :costS7]
        n = Symbol.(["QALYS$S" for S=1:7]) #[:QALYS1... :QALYS7]
        m = Symbol.(["Qclin$S" for S=1:7])
        R = Symbol.(["QGW$S" for S=1:7])
        P = Symbol.(["QICU$S" for S=1:7])
        append!(names,n)
        append!(names,m)
        append!(names,R)
        append!(names,P)
        df_table = DataFrame(table,names)
        CSV.write("CostEff_datatable_$season.csv",df_table)
    end
end
export CostEff_datatable


#-----------------
# get clean infection for calculation.
function get_inf(season::Symbol)
    #---------------------------
    # Reading total infection (results) from Agent-based Simulation
    # introduced_inf | inf0to2_total | inf0to2_H | inf0to2_P | ...500×14 DataFrame table
    infection = read_file("infection_$season.csv") # infections per simulation in a row
    _inf = select(infection, Not([1,2,5,8,11])) # take out 'total infection' columns per agegroup
    inf = Matrix(_inf) # convert to matrix

    ######### due to lack of hospitalization data for 2 and 3 years old, put zeros on their respective columns
    nrows, ncolm = size(inf)
    inf[:,[9,8,7]] = zeros(Int64,(nrows,3))

    return inf
end
export get_inf


#------------- calculated QALY for averted medically-atendent RSV infection
function Eff_averted(S::Int64,season::Symbol,inf::Array{Int64,2})

    # reading data file from my working directory
    _clinic = read_file("clinic_S$S$season.csv")
    clinic = Matrix(_clinic) # convert DataFrame to matrix

    #-------- reading data file from my working directory
    _GW = read_file("GW_S$S$season.csv")
    GW = Matrix(_GW) # convert DataFrame to Matrix

    _ICU = read_file("ICU_S$S$season.csv")
    ICU = Matrix(_ICU) # convert DataFrame to Matrix

    avert = inf - (clinic + GW + ICU)


    Q_avert = similar(avert,Float64)
    nrows, ncolm = size(Q_avert)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            n    = avert[sim,age]

            #calculate healthy QALY
            # disutilities obtained from
            #ref: Roy LM. Deriving health utility weights for infants with Respiratory Syncytial Virus (RSV)
            #(T). Univer- sity of British Columbia. Available: https://open.library.ubc.ca/cIRcle/collections/24/items/1.0074259.
            U2 = rand(Beta(8.41,1.15),n) # quality of life < 16 years old mean 0.88 HEALTHY
            Q_avert[sim,age] = sum(U2) # hospitalization
        end
    end

    # add up all age groups
    Q_avert_to = [sum(Q_avert[sim,:]) for sim=1:nrows]

 return Q_avert_to

end
export Eff_averted
#---------------------------------------------------
# Calculate total costs/days per simulation level
#---------------------------------------------------
function CostEff_total(S::Int64,season::Symbol)
    # get costs and effects for all agegroups and health status at the birth
    C_age , E_age, Q_clinic, Q_GW, Q_ICU = CostEff(S,season)
    nrows, ncolm = size(C_age)
    C_total = [sum(C_age[sim,:]) for sim=1:nrows]
    E_total = [sum(E_age[sim,:]) for sim=1:nrows]
    Q_clin  = [sum(Q_clinic[sim,:]) for sim=1:nrows]
    Q_W     = [sum(Q_GW[sim,:]) for sim=1:nrows]
    G_CU    = [sum(Q_ICU[sim,:]) for sim=1:nrows]

    return C_total, E_total, Q_clin, Q_W , G_CU
end
export CostEff_total



#---------------------------------------------------------------------------
# Calculate total costs and total hospital days for different agegroups
#---------------------------------------------------------------------------
function CostEff(S::Int64,season::Symbol)
    # calculate total costs per agegroup
    clinic, Q_clinic = cost_clinic(S,season)
    GW, ICU, Q_GW, Q_ICU = cost_HospICU(S,season)

    if S == 1
        cost_total = clinic + GW + ICU
    else
        prophylaxis = cost_palivi_LAMA_ResVax(S)
        cost_total = clinic + GW +  ICU + prophylaxis
    end


    # calculate total QALY from 1 year period
    Q_total = Q_clinic + Q_GW + Q_ICU

    # results
    return cost_total, Q_total,Q_clinic,Q_GW, Q_ICU
end
export CostEff


#----------------------------------------------------
# Calculate costs associate with patients in
# 1) clinic, 2)regional hosp + tertiary general ward, 3)tertiary ICU

#------------ Local Clinic ----------------------
function cost_clinic(S::Int64,season::Symbol)
    # average cost per case for mild RSV infection (no hospital admission)
    cost_mild = 1569

    # reading data file from my working directory
    _clinic = read_file("clinic_S$S$season.csv")
    clinic = Matrix(_clinic) # convert DataFrame to matrix

    inf_day = similar(clinic,Float64)
    Q_clinic = similar(clinic,Float64)


    # calculate total costs
    cost = clinic .* cost_mild

    # calculate infection days
    nrows, ncolm = size(clinic)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            n = clinic[sim,age] # number of patients
            day = get_inf_day(n) # total number of infection days per agegroup

            #calculate QALY
            # disutilities obtained from
            #ref: Roy LM. Deriving health utility weights for infants with Respiratory Syncytial Virus (RSV)
            #(T). University of British Columbia. Available: https://open.library.ubc.ca/cIRcle/collections/24/items/1.0074259.
            # inout is (day,mean,standard error)
            Q_clinic[sim,age] = sum(get_QALY(n,day,0.16,0.02)) # outpatient visits
        end
    end

    return cost, Q_clinic
end
export cost_clinic



#------- calculate infection days ---------------
function get_inf_day(n::Int64)
    expo = rand(Uniform(4.54,5.37),n) # exposure days
    symp = rand(Uniform(5.68,6.63),n) # symptomatic days
    inf_day = expo + symp

    return inf_day
end
export get_inf_day
#------------------ Regional Hospital ---------------
function cost_HospICU(S::Int64, season::Symbol)
    #------- Parameters
    # ------- AVERAGE COSTS per inpatient case ended up to hospital or ICU
    cost_GW = 16946
    cost_ICU = 80590

    # range for minimum stay in general ward or ICU per agegroups
    #GWStay = @SVector[1:8,1:8,1:7,1:7,1:10,1:10,  0:0,0:0,0:0]
    #ICUstay = @SVector[2:25,2:25,2:25,2:25,4:15,4:15,   0:0,0:0,0:0]

    GWStay = @SVector[3:3,3:3,3:3,3:3,4:4,4:4,  0:0,0:0,0:0]
    ICUstay = @SVector[12.5:12.5,12.5:12.5,12.5:12.5,12.5:12.5,9.5:9.5,9.5:9.5,   0:0,0:0,0:0]

    #-------- reading data file from my working directory
    _GW = read_file("GW_S$S$season.csv")
    GW = Matrix(_GW) # convert DataFrame to Matrix

    _ICU = read_file("ICU_S$S$season.csv")
    ICU = Matrix(_ICU) # convert DataFrame to Matrix



    #--------------------- total length of hospital days per agegroup
    # hospital and ICU stays per individual are selected and added up for each agegroup
    daysGW  = similar(GW ,Float64) #collecting array
    daysICU = similar(ICU,Float64)
    Q_GW    = similar(ICU,Float64)
    Q_ICU   = similar(ICU,Float64)

    nrows, ncolm = size(GW)
    @inbounds for sim=1:nrows
        Random.seed!(sim)
        for age=1:ncolm
            n_GW   = GW[sim,age]
            n_ICU  = ICU[sim,age]
            daysGW  = rand(GWStay[age],n_GW)
            daysICU = rand(ICUstay[age],n_ICU)

            #calculate QALY
            # disutilities obtained from
            #ref: Roy LM. Deriving health utility weights for infants with Respiratory Syncytial Virus (RSV)
            #(T). Univer- sity of British Columbia. Available: https://open.library.ubc.ca/cIRcle/collections/24/items/1.0074259.
            # inout is (day,mean,standard error)
            Q_GW[sim,age] = sum(get_QALY(n_GW,daysGW,0.41,0.03)) # hospitalization
            Q_ICU[sim,age]  = sum(get_QALY(n_ICU,daysICU,0.60,0.03)) # ICU admission
        end
    end



    #-------- total cost for all hospital/ICU admissions per agegroup
    costsGW = similar(GW,Float64) #collecting array
    costsGW = GW .* cost_GW

    costsICU = similar(ICU,Float64) #collecting array
    costsICU = ICU .* cost_ICU # sever infected cases

    return costsGW, costsICU, Q_GW, Q_ICU
end
export cost_HospICU


#---------calculate QALY
# mu is mean and se is standard error
function get_QALY(n,day,mu,se)
    alpha, beta = get_beta_parm(mu, se)
    disutil = rand(Beta(alpha,beta),n) # disutility
    U1 = 1 .- disutil # utility while affecting by RSV
    U2 = rand(Beta(8.41,1.15),n) # quality of life < 16 years old mean 0.88
    L = day ./ 365
    Q = U1 .* L + U2 .* (1 .- L)
    #sr = 0.015
    #Q = U1 .* (1+r).^L + U2 .* (1+r).^(1 .- L)
    #Q = U1 .* ((1 .- exp.(-r .* L)) ./ r)

    return Q
end
export get_QALY

# ------- estimate parameters for beta distribution
function get_beta_parm(mu::Float64, se::Float64)
    var = se^2
    alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta = alpha * (1 / mu - 1)
    return alpha, beta
end
export get_beta_parm

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
    if S ==4 || S == 7
        #ResVax only
        @inbounds for sim=1:nrow
            Random.seed!(sim)
            total[sim,1] = 1560 # healthy term < 3 months
            total[sim,2] = 1560 # preterm < 3 months
        end

        # ResVax + LAMA
        if S == 7
            # one doses for only preterm_ill > 3 months
            @inbounds for sim=1:nrow
                Random.seed!(sim)
                for j=2:l_p
                    total[sim,j+j] = paliviCost_p[j] #LAMA price same as paliv.
                end
            end
        end

    # palivi. or LAMA for high risk
    elseif S == 2 || S == 3
        @inbounds for j=1:l_p
            for sim=1:nrow
                Random.seed!(sim)
                # five doses for only preterm_ill
                # no cost for healthy terms; remained zero
                if S == 2 #Palivi
                    total[sim,j+j] = 5 .* paliviCost_p[j]
                elseif S == 3 #LAMA
                    total[sim,j+j] = paliviCost_p[j]
                end
            end
        end

    # palivi. or LAMA for high risk + healthy
    elseif S == 5 || S == 6
        @inbounds for sim=1:nrow
            Random.seed!(sim)
            if S == 5 #Palivi
                # cost for healthy terms < 3 months old
                total[sim,1] = 5 .* paliviCost_H
                for j=1:l_p
                    # five doses for only preterm_ill
                    total[sim,j+j] = 5 .* paliviCost_p[j]
                end
            elseif S == 6 #LAMA
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




#---------- NOT USED FUNCTIONS FOR RESULT

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
        CI_min = zeros(Float64,(l_0,l_1))
        CI_max = zeros(Float64,(l_0,l_1))
        @inbounds for S_0 in senario_0
            for S_1 = S_0+1:l_1
                if S_0==2 && S_1==6 continue end
                if S_0==3 && S_1==7 continue end

                #ICER with different bases plus CI intervals
                ICER_mean, ICER_CIlow, ICER_CIup = get_ICER(S_0,S_1,season)

                #collect data
                table[S_0,S_1] = ICER_mean #mean
                CI_min[S_0,S_1] = ICER_CIlow # 95% CI_min
                CI_max[S_0,S_1] = ICER_CIup # 95% CI_max
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
        CSV.write("ICER_mean_$season.csv",df_ICER_table)
        CSV.write("ICER_CImin_$season.csv",df_ICER_CImin)
        CSV.write("ICER_CImax_$season.csv",df_ICER_CImax)
    end #season
end
export main_ICER


#----------------------------------------------
# self-created bootstrapped and CI code (NOT USED FOR RESULTS)
function get_ICER(S_0::Int64,S_1::Int64,season::Symbol)
    #-------------- Total ICER per simulation
    # calculate bootstrapped incrementals
    inc_cost,inc_day = bootmean_incremental(S_0,S_1,season)

    # calculate ICERs per bootstrap mean of costs and days
    ICER_mean = inc_cost ./ inc_day

    m, CIlow, CIup = conf(ICER_mean; conf_level=0.95)

    """
    df_data = DataFrame(meanICERs=ICER_mean)
    CSV.write("ICER_S$S_0$S_1$season.csv",df_data)
    """
    return m, CIlow, CIup
end
export get_ICER

#------ CI
function conf(x; conf_level=0.95)
    # CI:95% - confidence interval
    l = length(x)
    m = mean(x)
    boot_delta = [m - x[i] for i=1:l]
    s = sort(boot_delta)
    c = conf_level
    delta_head = s[ceil(Int16,l * (1+c)/2)]
    delta_tail = s[floor(Int16,l * (1-c)/2)]

    CI_low = m - delta_head
    CI_high = m - delta_tail

    return m, CI_low, CI_high
end
export conf


# save bootstrapped data from cost-effectiveness Analysis
# 1- deltacosts 2-deltadays 3-ICER (ratio of costs/days) per scenario and per season
function main_bootmeanData()
    sea = [:mild,:moderate,:severe]  #level of outbreak depends on the season
    for season in sea
        senario_0 = @SVector[1,2,3]
        senario_1 = @SVector[1,2,3,4,5,6,7]
        l_0 = length(senario_0)
        l_1 = length(senario_1)

        #collecting file
        inc_cost = zeros(2000,15)
        inc_day = similar(inc_cost)
        i = 1
        @inbounds for S_0 in senario_0
            for S_1 = S_0+1:l_1
                inc_cost[:,i],inc_day[:,i] = bootmean_incremental(S_0,S_1,season)
                i += 1
            end # baseline
        end #alternative scenario

        names = Symbol.(["S0S1","S0S2","S0S3","S0S4","S0S5","S0S6",
        "S1S2","S1S3","S1S4","S1S5","S1S6",
        "S2S3","S2S4","S2S5","S2S6"])
        df_cost = DataFrame(inc_cost,names)
        df_day = DataFrame(inc_day,names)
        CSV.write("bootstrapped_deltacost_$season.csv",df_cost)
        CSV.write("bootstrapped_deltaday_$season.csv",df_day)
    end #season
end
export main_bootmeanData


function bootmean_incremental(S_0::Int64,S_1::Int64,season::Symbol)
    #-------------- Total incremental per simulation
    C1_total, E1_total = CostEff_total(S_1,season) #costs and effects for S1
    C0_total, E0_total = CostEff_total(S_0,season) #costs and effects for base (S0)

    # bootstrapping mean distribution by 2000 replications
    mean_cost1 = bootstrap_mean(C1_total,2000)
    mean_day1 = bootstrap_mean(E1_total,2000)
    mean_cost0 = bootstrap_mean(C0_total,2000)
    mean_day0 = bootstrap_mean(E0_total,2000)

    # calculate incrementals per bootstrap mean of costs and days
    inc_cost = mean_cost1 - mean_cost0
    #averting hospital day is a desired effect. therefore reduced days shall convert to positive
    inc_day = -(mean_day1 - mean_day0)


    #save data
    #df_data = DataFrame(delta_costs=inc_cost, delta_days=inc_day, ICER=ratio)
    #CSV.write("bootstrappedData_S$S_0$S_1$season.csv",df_data)
    return inc_cost, inc_day
end
export bootmean_incremental

"""
#----------------------------------------------------
# Covariance Confidence Surface estimator data pair (x,y)
#---------------------------------------------------
function ErrorEllipse(x::Array{Float64,1},y::Array{Float64,1})
    # calculate eigenvalues and eigenvectors
    covariance = [cov(x,x,corrected=true) cov(x,y,corrected=true);
                  cov(y,x,corrected=true) cov(y,y,corrected=true)]
    val = eigvals(covariance)
    vec = eigvecs(covariance)


    a=sqrt(5.991*maximum(val)) #major axis
    b=sqrt(5.991*minimum(val)) #minor axis
    # the ellipse in x and y coordinates
    theta_grid = range(0,stop=2*pi,length=10)
    ellipse_x_r  = a*cos(theta_grid)
    ellipse_y_r  = b*sin(theta_grid)
    y = a*cos(theta_grid) + b*sin(theta_grid)

end
export ErrorEllipse


clear all;
close all;

% Create some random data
s = [2 2];
x = randn(334,1);
y1 = normrnd(s(1).*x,1);
y2 = normrnd(s(2).*x,1);
data = [y1 y2];

% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
hold on;

% Plot the original data
plot(data(:,1), data(:,2), '.');
mindata = min(min(data));
maxdata = max(max(data));
Xlim([mindata-3, maxdata+3]);
Ylim([mindata-3, maxdata+3]);
hold on;

% Plot the eigenvectors
quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');

"""
#---------------------------------------------------------------------
# Bootstrap mean function for simulation outputs (costs,days)
#---------------------------------------------------------------------
function bootstrap_mean(A,n_boot::Int64)
    leng = length(A) # shall be number of simulations in row
    boot_mean = zeros(Float64,n_boot) # collecting array
    @inbounds for i=1:n_boot
        Random.seed!(i)
        # get the mean of a sample with replecementmean from our original data A
        m = mean(sample(A,leng; replace=true))
        boot_mean[i] = m # bootstrap mean distribution
    end
    return boot_mean
end
export bootstrap_mean

#------------------------------------------------------------------
# interval% CI - confidence interval accelerated bias-corrected method
#------------------------------------------------------------------
function confinterval(A,n_boot,interval::Float64) #interval shall be from [0,1]
    b = bootstrap(mean,A, BasicSampling(n_boot))
    c = confint(b,BCaConfInt(interval))
    m = c[1][1] # mean
    CI_lower = c[1][2] # lower bound of CI
    CI_upper = c[1][3] # upper bound of CI

    return m, CI_lower, CI_upper
end
export confinterval



#---------------------------------------------


#--------------------------------------
# to run in REPL
@everywhere using .RSV_CEA
####################
end #module
