## single agent-based model for Respiratory syncytial virus (RSV)
## Developed by Shokoofeh & Affan

module RSV_Results

## Loading necessary packages

using Distributed
#using Distributions
#using StatsBase
#using StaticArrays
using Random
using DataFrames
using CSV
#using Plots
#using StatsPlots

#using PyCall
#plt = pyimport("matplotlib.pyplot")
#sns = pyimport("seaborn")

# infection and hospitalized cases
# High/Medium/Mild outbreaks

function Boplot()
    #load prepared data
    dinf_high,dinf_mid,dinf_low, dhosp_high,dhosp_mid,dhosp_low=data_to_plot()
    #transform dataset from wide to long form
    infData = wide_to_long_data(dinf_high,dinf_mid,dinf_low)
    hospData = wide_to_long_data(dhosp_high,dhosp_mid,dhosp_low)

    #grouped boxplot
    #ax = sns.boxplot(x="age",y="value",data=infData)
end
export Boplot



function wide_to_long_data(d_high,d_mid,d_low)
    row,colm = size(d_high)
    age = ["0-2" "0-2" "3-5" "3-5" "6-11" "6-11"]
    Id = ["healthy-term" "preterm/ill" "healthy-term" "preterm/ill" "healthy-term" "preterm/ill"]
    ddd = [d_high,d_mid,d_low]
    level = ["high" "mid" "low"]

    # collecting table
    dataset = DataFrame([String,String,String,Float64],[:age,:healthId,:seasonLevel,:value])
    for j=1:3
        data = ddd[j]
        for i=1:colm
            d = DataFrame(age=age[i],healthId=Id[i],seasonLevel=level[j],value=data[:,i])
            append!(dataset,d)
        end
    end
    return dataset
end
export wide_to_long_data

function data_to_plot()
    #load simulation results
    infection = read_file("infection.csv")
    _inf = select(infection, Not([1,2,5,8,11,12,13,14]))
    hospital = read_file("HospS3.csv")
    _hosp = select(hospital, Not([7,8,9]))
    population = read_file("population.csv")
    _pop = select(population, Not([1,4,7,10,11,12,13]))

    # determine which simulations were high/low/mid season
    high = findall(i-> i.>=712 ,infection[:,1])
    mid = findall(i-> i.>550 && i.<712 ,infection[:,1])
    low = findall(i-> i.<=550 ,infection[:,1])

    # ------------------
    # infection
    inf_high = _inf[high,:]
    pop_high = _pop[high,:]
    hosp_high = _hosp[high,:]
    inf_mid = _inf[mid,:]
    pop_mid = _pop[mid,:]
    hosp_mid = _hosp[mid,:]
    inf_low = _inf[low,:]
    pop_low = _pop[low,:]
    hosp_low = _hosp[low,:]

    inf1000_high = (Matrix(inf_high)./Matrix(pop_high))*1000
    replace!(inf1000_high, NaN=>0)
    inf1000_mid = (Matrix(inf_mid)./Matrix(pop_mid))*1000
    replace!(inf1000_mid, NaN=>0)
    inf1000_low = (Matrix(inf_low)./Matrix(pop_low))*1000
    replace!(inf1000_low, NaN=>0)

    hosp1000_high = (Matrix(hosp_high)./Matrix(pop_high))*1000
    replace!(hosp1000_high, NaN=>0)
    hosp1000_mid = (Matrix(hosp_mid)./Matrix(pop_mid))*1000
    replace!(hosp1000_mid, NaN=>0)
    hosp1000_low = (Matrix(hosp_low)./Matrix(pop_low))*1000
    replace!(hosp1000_low, NaN=>0)

    return inf1000_high,inf1000_mid,inf1000_low,hosp1000_high,hosp1000_mid,hosp1000_low
end
export data_to_plot




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

"""
    boxplot(["0-2 months" "3-5 months" "6-11 months"],
            [ratio[:,2],ratio[:,4],ratio[:,6]],
            leg = true,
            colour= [:green :green :green],
            ylabel = "Infection (percentage)")


    plot!(legend=false, xtickfont=font(12), ytickfont=font(12),guidefont = font(15))
"""
#--------------------------------------
# to run in REPL
@everywhere using .RSV_Results
####################
end #module
