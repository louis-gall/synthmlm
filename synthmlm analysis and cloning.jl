using GLMakie, Flux, DataFrames, JLD2, CSV, FHist, Statistics, Distributions, LinearAlgebra, Integrals, ProgressMeter, SpecialFunctions, StaticArrays, Random

include("synthmlm base.jl")

#DATASET DOWNLOAD LINK:
#https://bham-my.sharepoint.com/personal/l_gall_bham_ac_uk/_layouts/15/guestaccess.aspx?share=IQCYxNtLeDAQS72oKgQc4K-CAapvFvztziluSq1UcCLLt8E&e=ncLM9o

#MODEL DOWNLOAD LINK:
#https://bham-my.sharepoint.com/personal/l_gall_bham_ac_uk/_layouts/15/guestaccess.aspx?share=IQDfzZ0WgJl1TquMVPdj3Hx3AT_b4Y-7QrOfMexzT-rHilI&e=D3bYmm

println("Make sure the model and the trial dataset are downloaded (details in readme and above this text)")
println("To use this file, you will need curvheuristics.jld in the same folder as this script.")
println("Thanks for downloading SynthMLM. This script demonstrates the dataset analysis and cloning component of SynthMLM")
println("----------")
println("Run analysedataset(location) to analyse a dataset stored at location.")
println("Run clonedataset(j,location) to clone a dataset stored at location, creating j synthetic images.")
println("Default parameters have the dataset location as nocodazole_1 (also in this Github repo), which when")
println("placed in the same folder as this file will allow the code to analyse and/or clone this data.")
println("You can analyse/clone your own data by replacing this folder. It will need to be formated in the way nano-org formats its data.")
println("----------")
println("Analysing and producing descriptor-matched synthetic clones of datasets is far easier on nano-org.bham.ac.uk")
println("----------")


function protomod(input,model,pmax,pmin)
    output = relu(model(input))

    for i in axes(output,2)
        for j in axes(output,1)
            output[j,i] = (pmax[j] - pmin[j])*output[j,i] + pmin[j]
        end
    end

    isdiscrete = [1 14 19 26]
    for i in isdiscrete
        output[i,:] = round.(Int,output[i,:])
    end

    for i in axes(output,2)
        if output[1,i] == .0
            output[2:13,i] .= 0
        end
        if output[14,i] == .0
            output[15:18,i] .= 0
        end
        if output[19,i] == .0
            output[20:25,i] .= 0
        end
        if output[26,i] == .0
            output[27:30,i] .= 0
        end
    end

    return output
end

function setupmodel()
    @load "synthmlm_analysis.jld2" model model_state pmax pmin
    Flux.loadmodel!(model, model_state)
    testmode!(model)

    modelpred(input) = protomod(input,model,pmax,pmin)

    return modelpred, length(pmax)
end

## END EXPLICITLY NEW STUFF

function csv2image(files)
    data = Array{Float32,4}(undef,100,100,1,length(files))
    for k in axes(data,4)
        temp = CSV.read(files[k],DataFrame)
        if in("x [nm]", names(temp))
            rename!(temp, Dict("x [nm]" => "x[nm]", "y [nm]" => "y[nm]"))
        end
        select!(temp,["x[nm]", "y[nm]"])
        temp = Matrix(temp)
        data[:,:,1,k] = bincounts(Hist2D((temp[:,1],temp[:,2]), binedges=(range(minimum(temp[:,1]), maximum(temp[:,1]), length=101),range(minimum(temp[:,2]), maximum(temp[:,2]), length=101))))
    end
    
    scaley = 40 
    clamp!(data,0,scaley)
    data ./= scaley

    return data
end

# example function that does a single dataset
function getanalysis(DIRECTORYNAME)
    modelpred, numvar = setupmodel()
    files = readdir(DIRECTORYNAME,join=true)
    filter!(x -> !occursin(".DS",x),files)
    filter!(x -> occursin(".csv",x),files)

    data = csv2image(files) # NAME CHANGE
    realdata = deepcopy(data)
    output = modelpred(data)
   
    return output, numvar, realdata
end

function getmeansandvars(i,output,hasfibs,hasclusts,hasrings,hasback)
    # ["numfibres","dens","denstd","width","widthstd","curv","curvstd","len","lenstd","delta","deltastd","sig","sigstd","numnuc_c","dupenum_c","dupestd_c","sigma_c","sigmastd_c", "numrings","radii","radstd","width_r","widthstd_r","dens_r","densstd_r"," numnuc_b","dupenum_b","dupestd_b","sigma_b","sigstd_b"]
    if i in 2:13
        if all(output[i,hasfibs] .== .0)
            return 0
        else
            return Float64(mean(output[i,hasfibs]))
        end
    elseif i in 15:18
        if all(output[i,hasclusts] .== .0)
            return 0
        else
            return Float64(mean(output[i,hasclusts]))
        end
    elseif i in 20:25
        if all(output[i,hasrings] .== .0)
            return 0
        else
            return Float64(mean(output[i,hasrings]))
        end
    elseif i in 27:30
        if all(output[i,hasback] .== .0)
            return 0
        else
            return Float64(mean(output[i,hasback]))
        end
    else
        if all(output[i,:] .== 0)
            return 0
        elseif count(output[i,:] .== 0)/length(output[i,:]) >= 0.9
            return 0
        else
            return round.(Int,output[i,:])
        end
    end
end

function getPDvalues(output,numvar)
    hasfibs = findall(output[1,:] .!= 0)
    hasclusts = findall(output[14,:] .!= 0)
    hasrings = findall(output[19,:] .!= 0)
    hasback = findall(output[26,:] .!= 0)

    PDvals = Vector{Any}(undef,numvar)

    for i in 1:numvar
        PDvals[i] = getmeansandvars(i,output,hasfibs,hasclusts,hasrings,hasback)
    end
    
    return PDvals
end

function getdists(PDvals)
    discrete = [1 14 19 26]
    means = [1 2 4 6 8 10 12 14 15 17 19 20 22 24 26 27 29]
    stds = [1 3 5 7 9 11 13 14 16 18 19 21 23 25 26 28 30]

    maxdupe = 8*PDvals[17]
    dupe_c_up = minimum([maxdupe, 1000])

    lowers = [0, 0.25, 5, 0, 500, 0.01*PDvals[8],0.5*PDvals[12],0,11,5,0,25,3,0.1,0,1,5]
    if 3*PDvals[10] < 0.5*PDvals[12] 
        sigup = 0.5*PDvals[12] 
    else
        sigup = 3*PDvals[10] 
    end

    uppers = [0, 6, 300, 20, 10000, 0.5*PDvals[8], sigup, 0, dupe_c_up, 300, 0, 300, Inf, 4.5, 0, 10, 180]
    
    for i in eachindex(lowers)
        if in(means[i], discrete)
        else
            if lowers[i] > PDvals[means[i]]
                lowers[i] = 0.9*PDvals[means[i]]
            end
        end
    end

    for i in eachindex(lowers)
        if PDvals[stds[i]] == .0
            lowers[i] = -Inf
        end
    end

    for i in eachindex(uppers)
        if uppers[i] < lowers[i]
            uppers[i] = lowers[i]
        end
    end

    dists = Vector{Any}(undef,17)

    if length(PDvals[1]) > 1
        dists[1] = truncated(MixtureModel(map(u->Poisson(u), PDvals[1])),lower = minimum(PDvals[1]), upper = maximum(PDvals[1]))
    else
        dists[1] = Normal(PDvals[1],PDvals[1]/4)
    end
    dists[2] = truncated(Normal(PDvals[2],PDvals[3]),lower = lowers[2], upper = uppers[2]) # dens
    dists[3] = truncated(Normal(PDvals[4],PDvals[5]),lower = lowers[3], upper= uppers[3]) # width
    dists[4] = truncated(Normal(PDvals[6],PDvals[7]),lower = lowers[4], upper= uppers[4]) # curv
    dists[5] = truncated(Normal(PDvals[8],PDvals[9]),lower = lowers[5], upper = uppers[5]) # len
    dists[6] = truncated(Normal(PDvals[10],PDvals[11]),lower = lowers[6],upper = uppers[6]) # del
    dists[7] = truncated(Normal(PDvals[12],PDvals[13]),lower = lowers[7] ,upper = uppers[7]) # sig

    if length(PDvals[14]) > 1
        dists[8] = PDvals[14] # truncated(MixtureModel(map(u->Poisson(u), PDvals[14])),lower = minimum(PDvals[14]), upper = maximum(PDvals[14]))
    else
        dists[8] = Normal(PDvals[14],PDvals[14]/4)
    end
    dists[9] = truncated(Normal(PDvals[15],PDvals[16]),lower=lowers[9],upper = uppers[9]) # dupe_c
    dists[10] = truncated(Normal(PDvals[17],PDvals[18]),lower=lowers[10],upper = uppers[10]) # sigma_c

    if length(PDvals[19]) > 1
        dists[11] = truncated(MixtureModel(map(u->Poisson(u), PDvals[19])),lower = minimum(PDvals[19]), upper = maximum(PDvals[19]))
    else
        dists[11] = Normal(PDvals[19],PDvals[19]/4)
    end
    dists[12] = truncated(Normal(PDvals[20],PDvals[21]),lower=lowers[12],upper = uppers[12])
    dists[13] = truncated(Normal(PDvals[22],PDvals[23]),lower=lowers[13])
    dists[14] = truncated(Normal(PDvals[24],PDvals[25]),lower=lowers[14],upper = uppers[14])
    
    if length(PDvals[26]) > 1
        dists[15] = PDvals[26] # truncated(MixtureModel(map(u->Poisson(u), PDvals[26])),lower = minimum(PDvals[26]), upper = maximum(PDvals[26]))
    else
        dists[15] = Normal(PDvals[26],PDvals[26]/4)
    end
    dists[16] = truncated(Normal(PDvals[27],PDvals[28] ),lower=lowers[16],upper = uppers[16]) # dupe_b
    dists[17] = truncated(Normal(PDvals[29],PDvals[30] ),lower=lowers[17],upper = uppers[17]) # sigma_b
    
    return dists
end

function getparameterdistributions(output,numvar)
    # relationship between measured length and initial length in algorithm
    output[7,:] .*= 2
    output[8,:] .*= 2
    
    PDvals = getPDvalues(output,numvar)

    dists = getdists(PDvals)

    return dists
end

function getbenddeg(curv) # this has a file, basically I have made an inverse function from data I collected from the model output.
    @load "curvheuristics.jld2" curvs bends bezdegs
    # just limit degree?
    if curv < 15
        indies = findall(x-> x >= curv - 0.5 && x < curv + 0.5, curvs)
        bends = rand(bends[indies])
        bezdegs = maximum([rand(bezdegs[indies]), 2])
    else
        indies = findall(x-> x > 15, curvs)
        bends = rand(bends[indies])
        bezdegs = maximum([rand(bezdegs[indies]), 2])
    end

    return bends, bezdegs
end

function generatesynthdata(dists,num=10,inp_name="synthtest")
    if num == 1
    elseif mod(num,10) != 0
        println("Rounding number to nearest 10")
        num = round(Int, num / 10) * 10
    end

    data = Array{Float32,4}(undef,100,100,1,num)
    for i in axes(data,4)
        if dists[1] == 0
            initn = 0
        else
            initn = rand(dists[1])
        end
        if initn == 0
            lendist = Bernoulli(0)
            densdist_f = Bernoulli(0)
            widthdist_f = Bernoulli(0)
            deldist = Bernoulli(0)
            sigdist = Bernoulli(0)
            bend = 0
            bezdeg = 0
        else
            lendist = dists[5] 
            densdist_f = dists[2]
            widthdist_f = dists[3]
            deldist = dists[6]
            sigdist = dists[7]
            curv = rand(dists[4])
            bend, bezdeg = getbenddeg(curv)
        end

        if dists[8] == 0
            numnuc_c = 0
        else
            numnuc_c = rand(dists[8])
        end
        if numnuc_c == 0
            dupedist_c = Bernoulli(0)
            sigmadist_c = Bernoulli(0)
        else
            dupedist_c = dists[9]
            sigmadist_c = dists[10]
        end

        if dists[11] == 0
            initnumrings = 0
        else
            initnumrings = rand(dists[11])
        end
        if initnumrings == 0
            raddist = Bernoulli(0)
            widthdist_r = Bernoulli(0)
            densdist_r = Bernoulli(0)
        else
            raddist = dists[12]
            widthdist_r = dists[13]
            densdist_r = dists[14]
        end

        if dists[15] == 0
            numnuc_b = 0
        else
            numnuc_b = rand(dists[15])
        end
        if numnuc_b == 0
            dupedist_b = Bernoulli(0)
            sigmadist_b = Bernoulli(0)
        else
            dupedist_b = dists[16]
            sigmadist_b = dists[17]
        end
        bezdeg = 2
        data[:,:,1,i], _, _, _, _, _, _, _, _, _, _, _, _, _ = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,true,inp_name)
     end
    
    if num > 1
        ROI = [0 0; 30 0; 30 3*ceil(Int,num/10); 0 3*ceil(Int,num/10)]
    else
        ROI = [0 0; 3 0; 3 3; 0 3]
    end

    CSV.write(string(inp_name,"_ROI_1.csv"), DataFrame(ROI,["X", "Y"]))
    return data
end

function plotimage(data,num,name) # this makes a 100x100 image that the model uses, and saves it as synthtest.png
    fig = Figure(size=(1000,800))
    axs = Array{Any,2}(undef,5,4)
    plotme = randperm(size(data,4))
    for i in 1:(minimum([20 num]))
        axs[i] = Axis(fig[fld1(i,5),mod1(i,5)])
        hidedecorations!(axs[i])
        hidespines!(axs[i])
        heatmap!(axs[i],data[:,:,1,plotme[i]],colorrange=[0,0.66])
    end

    display(GLMakie.Screen(),fig)
    save(string(name,"_image.png"),fig) # will eventually need a way of getting a generated name instead of just test
end

function analysedataset(location = "nocodazole_1")
    output, numvar, realdata = getanalysis(location)

    fig = Figure(size=(2000,1000),fontsize=32)
    titles = ["Number Fibres","Fibre Density","Fibre Width","Curvature","Fibre Length","Delta","Sigma", "Number Clusters","Particles Per Cluster", "Cluster Width", "Number Rings","Ring Radius","Ring Width","Ring Density", "Small Cluster Number","Particles per Small Cluster","Small Cluster Width"]
    
    plotme = [1 2 4 6 8 10 12 14 15 17 19 20 22 24 26 27 29]
    
    axs = Vector{Any}(undef,30)
    for i in 1:17
        axs[i] = Axis(fig[fld1(i,6),mod1(i,6)],title=titles[i],xticklabelsize=10,yticklabelsize=10,titlesize=14)
    end

    for j in 1:17
        hist!(axs[j],output[plotme[j],:],normalization=:pdf)
    end

    display(fig)
end

function clonedataset(j=20,location = "nocodazole_1")
    output, numvar, realdata = getanalysis(location)
    dists = getparameterdistributions(output,numvar)
    
    data = generatesynthdata(dists,j,string(location,"_clone"))
    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    plotimage(realdata,j,string(location,"_real"))
    plotimage(data,j,string(location,"_synth"))
end
;