## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh

module RSV_plots

## Loading necessary packages

using Distributed
#using Distributions
#using StatsBase
#using StaticArrays
using Random
using DataFrames
using CSV
using StatsPlots #import Plots too
using ColorBrewer
using Statistics
#using Plots


function total_caseICUGW()
    #load simulation results per senario and season
    season = [:mild, :moderate, :severe]
    for sea in season

        #collecting matrices for new dataset for plotting
        aveCases_ICU = zeros(Float64,1,7)
        aveCases_GW = similar(aveCases_ICU,Float64)
        #---- load hospital/ICU/GW results per scenario
        for S=1:7
            hospital = read_file("HospS$S$sea.csv")
            _hosp = select(hospital, Not([7,8,9]))
            hosp = Matrix(_hosp)

            intensiveCU = read_file("ICUS$S$sea.csv")
            _ICU = select(intensiveCU, Not([7,8,9]))
            ICU = Matrix(_ICU)

            GW = hosp - ICU #general ward

            # add up number per age group and then get average to obtain mean of total cases
            aveCases_GW[1,S] = mean(sum(GW,dims=2))
            aveCases_ICU[1,S] = mean(sum(ICU,dims=2))
        end #S

        #---------- save data
        # DataFrame to collect data
        # names
        names = Symbol.(["S0","S1","S2","S3","S4","S5","S6"])

        df_GW = DataFrame(aveCases_GW, names)
        df_ICU = DataFrame(aveCases_ICU,names)
        CSV.write("aveCase_GW_$sea.csv",df_GW)
        CSV.write("aveCase_ICU_$sea.csv",df_ICU)
    end #season

    return aveCases_GW, aveCases_ICU
end
export total_caseICUGW



function acceptableICER()
    costGW = 16946
    costICU = 80590
    season = [:mild, :moderate, :severe]

    #collecting arrays
    accep_ICER = zeros(Int64,3,3)
    for ss= 1:length(season)
        sea = season[ss]
        # ----------load data
        # average of total cases per all seven Scenarios
        # average of total days
        _GW = read_file("aveCase_GW_$sea.csv")
        GW = Matrix(_GW)
        _ICU = read_file("aveCase_ICU_$sea.csv")
        ICU = Matrix(_ICU)
        _days = read_file("days$sea.csv")
        days500 = Matrix(_days) #include 500 simularion
        days = mean(days500,dims=1) #average days

        # obtain threshold ICERs for base scenarios S0,S1 and S2
        threshold = (GW[1:3] .* costGW + ICU[1:3] .* costICU) ./ days[1:3]
        accep_ICER[ss,:] = [ceil(Int64,threshold[i]) for i=1:3]
    end

    return accep_ICER
end
export acceptableICER
# ---------------

function costEff_plane()
    season = [:mild, :moderate, :severe]
    #season = [:mild]
    for ss =1:length(season)
        sea = season[ss]
        # load infection results
        _delta_cost = read_file("bootstrapped_deltacost_$sea.csv")
        delta_cost = Matrix(_delta_cost)
        _delta_day = read_file("bootstrapped_deltaday_$sea.csv")
        delta_day = Matrix(_delta_day)

        #Average of incerementals
        ave_C = mean(delta_cost,dims=1)
        ave_DeltaCost = [ave_C[i] for i=1:length(ave_C)] #convert to vector
        ave_D = mean(delta_day,dims=1)
        ave_DeltaDay = [ave_D[j] for j=1:length(ave_D)] #convert to vector

        #Plotting
        col = []
        colS0 = [:blue for i=1:6]
        colS1 = [:green for i=1:5]
        colS2 = [:red for i=1:4]
        append!(col,colS0)
        append!(col,colS1)
        append!(col,colS2)

        # calculate acceptable ICER (willingness to pay)
        WTP = acceptableICER()

        h = ["S0S1","S0S2","S0S3","S0S4","S0S5","S0S6",
        "S1S2","S1S3","S1S4","S1S5","S1S6",
        "S2S3","S2S4","S2S5","S2S6"]
        x = ave_DeltaDay
        y = ave_DeltaCost ./10^5
        k = reshape(WTP[ss,:],(1,3))

        yline = (k .*x)./10^5

        # create p1,p2,p3
        @eval $(Symbol(:p,ss)) = plot(x,y,seriestype=:scatter,
        framestyle = :origin,
        series_annotations = text.(h, :bottom, 7, col),
        color=col)

        @eval $(Symbol(:p,ss)) = plot!(x,yline, legend=false, color=[:blue :green :red])
    end

    #layout
    plot(p1,p2,p3,layout = (1,3))
    """
    plot!(xlabel = "reduced hospital days",
    ylabel = "incremental costs",
    xtickfont=font(10), ytickfont=font(10),guidefont = font(13))
    col = [:blue,:blue,:blue,:blue,:blue,:blue,
    :green,:green,:green,:green,:green,
    :red,:red,:red,:red]
    """
end
export costEff_plane


function bar_plot()
    m_clinic,m_GW,m_ICU = clinicGWICU_to_plot()
    blu = ColorBrewer.palette("Blues",8) # manual color choices
    pur = ColorBrewer.palette("Purples",8)
    ora = ColorBrewer.palette("Oranges",7)
    red = ColorBrewer.palette("RdPu",7)

    groupedbar([m_GW[:,1] m_GW[:,3] m_GW[:,5] m_ICU[:,1] m_ICU[:,3] m_ICU[:,5]],
                bar_position =:dodge,
                bar_width=0.6,
                color = [blu[3] blu[5] blu[7] pur[3] pur[5] pur[7]],
                #color = [ora[2] ora[4] ora[5] red[2] red[4] red[5]],
                label=["0-2 GW" "3-5 GW" "6-11 GW" "0-2 ICU" "3-5 ICU" "6-11 ICU"],
                ylabel = "per 100 population")

    plot!(legend=true, xtickfont=font(1), ytickfont=font(12),
    guidefont = font(15))
end
export bar_plot

function infbarplot()
    inf_h = read_file("infDataToPlot_healthy.csv")
    inf_p = read_file("infDataToPlot_pretermill.csv")

    blu = ColorBrewer.palette("Blues",8) # manual color choices
    ora = ColorBrewer.palette("Oranges",7)
    #pur = ColorBrewer.palette("Purples",8)
    #red = ColorBrewer.palette("RdPu",7)

    groupedbar([mean(inf_h[1]) mean(inf_p[1]) mean(inf_h[4]) mean(inf_p[4]) mean(inf_h[7]) mean(inf_p[7])],
                bar_position =:dodge,
                bar_width=0.6,
                color = [blu[3] ora[2] blu[3] ora[2] blue[3] ora[2]],
                #color = [ora[2] ora[4] ora[5] red[2] red[4] red[5]],
                label=["Healthy" "Preterm/chronically-ill"],
                ylabel = "Infection per 100 population")


    plot!(legend=true, xtickfont=font(1), ytickfont=font(12),
    guidefont = font(15))
end
export infbarplot






function PWICUdata_to_plot()
    #collecting matrix per agegroup and health status
    m_clinic = zeros(7,6)
    m_GW = zeros(7,6)
    m_ICU = zeros(7,6)

    #load files
    population = read_file("population.csv")
    pop = select(population, Not([1,4,7,10,11,12,13]))

    for S=1:5
        _clinic = read_file("clinicS$S.csv")
        clinic = select(_clinic, Not([7,8,9]))
        _Hosp = read_file("HospS$S.csv")
        Hosp = select(_Hosp, Not([7,8,9]))
        _ICU = read_file("ICUS$S.csv")
        ICU = select(_ICU, Not([7,8,9]))
         GW = Matrix(Hosp) - Matrix(ICU)

        #cases per 100 population
        clinic100 = (Matrix(clinic)./Matrix(pop))*100
        replace!(clinic100, NaN=>0)
        GW100 = (Matrix(GW)./Matrix(pop))*100
        replace!(GW100, NaN=>0)
        ICU100 = (Matrix(ICU)./Matrix(pop))*100
        replace!(ICU100, NaN=>0)

        m_clinic[S,:] = median(clinic100,dims=1)
        m_GW[S,:] = median(GW100,dims=1)
        m_ICU[S,:] = median(ICU100,dims=1)
    end

    #LAMA intervention has same efficacy of palivizumab
    # S6 and S7 are the same ICU and clinics as S2 and S3
    m_clinic[6,:] = m_clinic[2,:]
    m_clinic[7,:] = m_clinic[3,:]
    m_GW[6,:] = m_GW[2,:]
    m_GW[7,:] = m_GW[3,:]
    m_ICU[6,:] = m_ICU[2,:]
    m_ICU[7,:] = m_ICU[3,:]

    return m_clinic,m_GW,m_ICU
end
export clinicGWICUdata_to_plot






"""
function dplot()
    inf_high,inf_mid,inf_low, hosp_high,hosp_mid,hosp_low=data_to_plot()
    for i=1:3
        indx = 2*i-1 #healthy term infants
        #indx = 2*i #preterm_ill
        StatsPlots.dotplot!([inf_high[:,indx],inf_mid[:,indx],inf_low[:,indx]],
                markersize=2,
                color = [blues[7] blues[5] blues[2]])
    end
end
export dplot
"""

function hospplots()
    hosp1,hosp2,hosp3,hosp4,hosp5 = hosp_to_plot()

    #blues = ColorBrewer.palette("Blues",7) # manual color choices
    oranges = ColorBrewer.palette("Oranges",7)
    for i=1:3 #three age category under 1 years old
        #indx = 2*i-1 #healthy term infants
        indx = 2*i #preterm_ill
        StatsPlots.boxplot!([hosp1[:,indx],hosp2[:,indx],hosp3[:,indx],hosp4[:,indx],hosp5[:,indx],hosp2[:,indx],hosp3[:,indx]],
                leg = true,
                outliers = true,
                #color = [blues[7] blues[6] blues[5] blues[4] blues[3] blues[2] blues[1]],
                color = [oranges[7] oranges[6] oranges[5] oranges[4] oranges[3] oranges[2] oranges[1]],
                bar_width = 0.4,
                ylabel = "hospital admissions per 100 popu.")
    end
    plot!(legend=false, xtickfont=font(1), ytickfont=font(12),guidefont = font(15))

end
export hospplots

function hosp_plots()
    # load infection results
    hosp_mild= read_file("hospDataToPlot_TOTAL_mild.csv")
    hosp_mod= read_file("hospDataToPlot_TOTAL_moderate.csv")
    hosp_sev= read_file("hospDataToPlot_TOTAL_severe.csv")

    #blues = ColorBrewer.palette("Blues",7) # manual color choices
    oranges = ColorBrewer.palette("Oranges",7)

    StatsPlots.boxplot!([hosp_sev[:,1] hosp_sev[:,2] hosp_sev[:,3] hosp_sev[:,4] hosp_sev[:,5] hosp_sev[:,6] hosp_sev[:,7]],
    leg = true,
    outliers = true,
    #color = [blues[7] blues[6] blues[5] blues[4] blues[3] blues[2] blues[1]],
    #color = [oranges[7] oranges[6] oranges[5] oranges[4] oranges[3] oranges[2] oranges[1]],
    bar_width = 0.4,
    ylabel = "hospital admissions per 100 pop.")

    plot!(legend=false, xtickfont=font(1), ytickfont=font(12),guidefont = font(15))

end
export hosp_plots

label=["S0" "S1" "S2" "S3" "S4" "S5" "S6"]


function hospICUGWdata_to_plot()
    #load simulation results per senario and season
    # AND then compute (hospital per infection)* 100
    season = [:mild, :moderate, :severe]
    for sea in season
        """
        # load infection results
        infection = read_file("infection_$sea.csv")
        _inf = select(infection, Not([1,2,5,8,11,12,13,14]))
        inf = Matrix(_inf)
        """
        #load population
        population = read_file("population.csv")
        _pop = select(population, Not([1,4,7,10,11,12,13]))
        pop = Matrix(_pop)

        #collecting matrices for new dataset for plotting
        hosp100_H = zeros(Float64,500,21)
        hosp100_P = similar(hosp100_H,Float64)
        ICU100_H = similar(hosp100_H,Float64)
        ICU100_P = similar(hosp100_H,Float64)
        GW100_H = similar(hosp100_H,Float64)
        GW100_P = similar(hosp100_H,Float64)
        #---- load hospital/ICU/GW results per scenario
        for S=1:7
            hospital = read_file("Hosp_S$S$sea.csv")
            _hosp = select(hospital, Not([7,8,9]))
            hosp = Matrix(_hosp)

            intensiveCU = read_file("ICU_S$S$sea.csv")
            _ICU = select(intensiveCU, Not([7,8,9]))
            ICU = Matrix(_ICU)

            GW = hosp - ICU #general ward

            #hospital/ICU/GW admissions per 100 infection
            hosp100 = (hosp./pop)*100
            ICU100 = (ICU./hosp)*100
            GW100 = (GW./hosp)*100

            """
            # remove NAN simulations where population is zero in denominator
            deleteat!(hosp100, hosp100.== NaN)
            deleteat!(ICU100, ICU100.== NaN)
            deleteat!(GW100, GW100.== NaN)
            """

            #create new Dataset for plotting (seperating H and P, combining intervention senarios (S) in one file)
            itr = [S,S+7,S+14] #index for new dataset
            @inbounds for col=1:3
                col_h = 2*col-1 #iteration of old data (H)
                col_p = 2*col #iteration of old data (P)
                c = itr[col]
                hosp100_H[:,c] = hosp100[:,col_h]
                hosp100_P[:,c] = hosp100[:,col_p]

                ICU100_H[:,c] = ICU100[:,col_h]
                ICU100_P[:,c] = ICU100[:,col_p]

                GW100_H[:,c] = GW100[:,col_h]
                GW100_P[:,c] = GW100[:,col_p]
            end
        end #S

        #---------- save data
        # DataFrame to collect data
        # names
        names = Symbol.(["0-2_S1","0-2_S2","0-2_S3","0-2_S4","0-2_S5","0-2_S6","0-2_S7",
                    "3-5_S1","3-5_S2","3-5_S3","3-5_S4","3-5_S5","3-5_S6","3-5_S7",
                    "6-11_S1","6-11_S2","6-11_S3","6-11_S4","6-11_S5","6-11_S6","6-11_S7"])

        df_hosp100H = DataFrame(hosp100_H, names)
        df_hosp100P = DataFrame(hosp100_P,names)
        df_ICU100H = DataFrame(ICU100_H, names)
        df_ICU100P = DataFrame(ICU100_P,names)
        df_GW100H = DataFrame(GW100_H, names)
        df_GW100P = DataFrame(GW100_P,names)
        CSV.write("hospDataToPlot_healthy_$sea.csv",df_hosp100H)
        CSV.write("hospDataToPlot_preterm_$sea.csv",df_hosp100P)
        """
        CSV.write("ICUDataToPlot_healthy_$sea.csv",df_ICU100H)
        CSV.write("ICUDataToPlot_preterm_$sea.csv",df_ICU100P)
        CSV.write("GWDataToPlot_healthy_$sea.csv",df_GW100H)
        CSV.write("GWDataToPlot_preterm_$sea.csv",df_GW100P)
        """
    end #season
end #function
export hospICUGWdata_to_plot

#-------------------------------------------------
# TOTAL means we add up data for all age categories under one group of 12 month old
function hospICUGWdata_to_plot_TOTAL()
    #load simulation results per senario and season
    # AND then compute (hospital per infection)* 100
    season = [:mild, :moderate, :severe]
    for sea in season
        # load infection results
        infection = read_file("infection_$sea.csv")
        _inf = select(infection, Not([1,2,5,8,11,12,13,14]))
        _inf_h = select(_inf,Not(2,4,6)) #healthy
        _inf_p = select(_inf,Not(1,3,5)) #preterm
        # add up age category
        infH = sum(Matrix(_inf_h),dims=2)
        infP = sum(Matrix(_inf_p),dims=2)

        #collecting matrices for new dataset for plotting
        hosp100_H = zeros(Float64,500,7)
        hosp100_P = similar(hosp100_H,Float64)
        ICU100_H = similar(hosp100_H,Float64)
        ICU100_P = similar(hosp100_H,Float64)
        GW100_H = similar(hosp100_H,Float64)
        GW100_P = similar(hosp100_H,Float64)
        #---- load hospital/ICU/GW results per scenario
        for S=1:7
            hospital = read_file("Hosp_S$S$sea.csv")
            _hosp = select(hospital, Not([7,8,9]))
            _hosp_h = select(_hosp,Not(2,4,6)) #healthy
            _hosp_p = select(_hosp,Not(1,3,5)) #preterm
            # add up all age category
            hospH = sum(Matrix(_hosp_h),dims=2)
            hospP = sum(Matrix(_hosp_p),dims=2)

            intensiveCU = read_file("ICU_S$S$sea.csv")
            _ICU = select(intensiveCU, Not([7,8,9]))
            _ICU_h = select(_ICU,Not(2,4,6))
            _ICU_p = select(_ICU,Not(1,3,5))
            # add up age category
            ICUH = sum(Matrix(_ICU_h),dims=2)
            ICUP = sum(Matrix(_ICU_p),dims=2)

            GWH = hospH - ICUH #general ward (healthy)
            GWP = hospP - ICUP #preterm

            #hospital admissions per 100 infection
            hosp100_h = (hospH./infH)*100
            hosp100_p = (hospP./infP)*100
            #replace!(hosp100_h, NaN=>0)
            #replace!(hosp100_p, NaN=>0)

            ICU100_h = (ICUH./hospH)*100
            ICU100_p = (ICUP./hospP)*100
            #replace!(ICU100_h, NaN=>0)
            #replace!(ICU100_p, NaN=>0)

            GW100_h = (GWH./hospH)*100
            GW100_p = (GWP./hospP)*100
            #replace!(GW100_h, NaN=>0)
            #replace!(GW100_p, NaN=>0)

            # keep data per loop of scenario S
            hosp100_H[:,S] = hosp100_h
            hosp100_P[:,S] = hosp100_p

            ICU100_H[:,S] = ICU100_h
            ICU100_P[:,S] = ICU100_p

            GW100_H[:,S] = GW100_h
            GW100_P[:,S] = GW100_p
        end #S

        #---------- save data
        # DataFrame to collect data
        # names
        names = Symbol.(["S0","S1","S2","S3","S4","S5","S6"])

        df_hosp100H = DataFrame(hosp100_H, names)
        df_hosp100P = DataFrame(hosp100_P,names)
        df_ICU100H = DataFrame(ICU100_H, names)
        df_ICU100P = DataFrame(ICU100_P,names)
        df_GW100H = DataFrame(GW100_H, names)
        df_GW100P = DataFrame(GW100_P,names)
        """
        CSV.write("hospDataToPlot_healthy_$sea.csv",df_hosp100H)
        CSV.write("hospDataToPlot_preterm_$sea.csv",df_hosp100P)
        """
        CSV.write("ICUDataToPlot_healthy_$sea.csv",df_ICU100H)
        CSV.write("ICUDataToPlot_preterm_$sea.csv",df_ICU100P)
        CSV.write("GWDataToPlot_healthy_$sea.csv",df_GW100H)
        CSV.write("GWDataToPlot_preterm_$sea.csv",df_GW100P)
    end #season
end #function
export hospICUGWdata_to_plot_TOTAL

#-------------------------------------------------

##################################################
# TOTAL means we add up data for all age categories under one group of 12 month old
function hospICUGWdata_to_plot_total()
    #load simulation results per senario and season
    # AND then compute hospital per 100 population

    #load population
    population = read_file("population.csv")
    _pop = select(population, [1,4,7])
    # add up age category
    pop = sum(Matrix(_pop), dims=2)

    ##
    season = [:mild, :moderate, :severe]
    for sea in season
        #collecting matrices for new dataset for plotting
        hosp100 = zeros(Float64,500,7)
        ICU100 = similar(hosp100,Float64)
        GW100 = similar(hosp100,Float64)

        #---- load hospital/ICU/GW results per scenario
        for S=1:7
            hospital = read_file("Hosp_S$S$sea.csv")
            _hosp = select(hospital, Not([7,8,9]))
            # add up all age category
            hosp = sum(Matrix(_hosp),dims=2)

            intensiveCU = read_file("ICU_S$S$sea.csv")
            _ICU = select(intensiveCU, Not([7,8,9]))
            # add up age category
            ICU = sum(Matrix(_ICU),dims=2)

            GW = hosp - ICU #general ward

            # calculate incidence and keeping them per loop of scenario S
            hosp100[:,S] = (hosp./pop)*100 #per 100 population
            ICU100[:,S] = (ICU./hosp)*100 #per 100 hospitalization
            GW100[:,S] = (GW./hosp)*100 #per 100 hospitalization
        end #S

        #---------- save data
        # DataFrame to collect data
        # names
        names = Symbol.(["S0","S1","S2","S3","S4","S5","S6"])

        df_hosp100 = DataFrame(hosp100, names)
        df_ICU100 = DataFrame(ICU100, names)
        df_GW100 = DataFrame(GW100, names)

        CSV.write("hospDataToPlot_TOTAL_$sea.csv",df_hosp100)
        CSV.write("ICUDataToPlot_TOTAL_$sea.csv",df_ICU100)
        CSV.write("GWDataToPlot_TOTAL_$sea.csv",df_GW100)
    end #season
end #function
export hospICUGWdata_to_plot_total


#save 'hosp100' in created files with names
#hosp100S(from 1 to 5)(mild or moderate or severe)
#@eval $(Symbol(:hosp100S,S)) = hosp100

function infplots(preterm::Bool)
    if preterm == true
        filename = "infDataToPlot_pretermill.csv"
        oranges = ColorBrewer.palette("Oranges",7)
        colorArray = [oranges[7] oranges[5] oranges[2]]
    else
        filename = "infDataToPlot_healthy.csv"
        blues = ColorBrewer.palette("Blues",7) # manual color choices
        colorArray = [blues[7] blues[5] blues[2]]
    end

    #load infection data for plotting, created by infData_to_plot()
    data = read_file(filename)

    itr = [1,4,7]
    for col in itr
        StatsPlots.boxplot!([data[:,col+2] data[:,col+1] data[:,col]],
                leg = true,
                outliers = true,
                color = colorArray,
                bar_width = 0.4,
                label=["severe season" "moderate season" "mild season"],
                #ylabel = "infection per 100 at-risk infants in age category"
                )
    end
    plot!(legend=false, xtickfont=font(1), ytickfont=font(12),guidefont = font(15))
end
export infplots

function infdata_to_plot()
    #load simulation results
    population = read_file("population.csv")
    _pop = select(population, Not([1,4,7,10,11,12,13]))
    pop = Matrix(_pop)

    infection_mild = read_file("infection_mild.csv")
    _infmild = select(infection_mild, Not([1,2,5,8,11,12,13,14]))
    infmild = Matrix(_infmild)

    infection_moderate = read_file("infection_moderate.csv")
    _infmod = select(infection_moderate, Not([1,2,5,8,11,12,13,14]))
    infmod = Matrix(_infmod)

    infection_severe = read_file("infection_severe.csv")
    _infsev = select(infection_severe, Not([1,2,5,8,11,12,13,14]))
    infsev = Matrix(_infsev)

    # ------------------
    #total infection per total population in every age category * 1000
    infmild100 = (infmild./pop)*100
    infmod100 = (infmod./pop)*100
    infsev100 = (infsev./pop)*100

    #create new Dataset for plotting (seperating H and P, combining seasons in one file)
    #collecting matrices
    infH = zeros(500,9)
    infP = similar(infH,Float64)

    # HEALTHY INFANTS
    itr_c_h =[1,4,7]
    itr_col_h =[1,3,5]
    @inbounds for i=1:length(itr_c_h)
        c = itr_c_h[i]
        col = itr_col_h[i]
        infH[:,c] = infmild100[:,col]
        infH[:,c+1] = infmod100[:,col]
        infH[:,c+2] = infsev100[:,col]
    end

    # PRETERM/ILL INFANTS
    itr_c_p =[1,4,7]
    itr_col_p =[2,4,6]
    @inbounds for i=1:length(itr_c_p)
        c = itr_c_p[i]
        col = itr_col_p[i]
        infP[:,c] = infmild100[:,col]
        infP[:,c+1] = infmod100[:,col]
        infP[:,c+2] = infsev100[:,col]
    end

    #---------- save data
    # DataFrame to collect data
    # names
    names = Symbol.(["0-2_mild","0-2_mod","0-2_severe",
                           "3-5_mild","3-5_mod","3-5_severe",
                           "6-11_mild","6-11_mod","6-11_severe"])

    df_infH = DataFrame(infH, names)
    df_infP = DataFrame(infP,names)
    CSV.write("infDataToPlot_healthy.csv",df_infH)
    CSV.write("infDataToPlot_pretermill.csv",df_infP)

    return
end
export infdata_to_plot




#----------------------------------------------
# reading csv.files from saved-data of the simulation
#-----------------------------------------------
function read_file(file::String)
    fname = file  #.csv file
    fpath = "/Users/shokoofehnourbakhsh/Dropbox/YorkUniversity/RSV_Shokoofeh/JuliaCodes/"
    #fpath = "/Users/shokoofehnourbakhsh/Dropbox/YorkUniversity/RSV_Shokoofeh/JuliaCodes/RSV_DatasetNOV2020/DataForPlot/"
    path = fpath*fname
    f = CSV.read(path,normalizenames=true) #DataFrame table
    return f
end
export read_file

#--------------------------------------
# to run in REPL
@everywhere using .RSV_plots
####################
end #module
