using GLMakie, Flux, DataFrames, JLD2, CSV, FHist, Statistics, Distributions, LinearAlgebra, Integrals, ProgressMeter, SpecialFunctions, StaticArrays, Random

include("synthmlm base.jl")

import Base.rand

# PRESS RUN ON THIS FILE TO COMPILE FUNCTIONS, THEN TYPE THE FUNCTIONS DESCRIBED BELOW IN TERMINAL
println("To use this file, you will need curvheuristics.jld in the same folder as this script.")
println("Thanks for downloading SynthMLM. This script has a few functions that demonstrate the sythetic data generation component of SynthMLM")
println("----------")
println("generateall(j,savepp) generates j synthetic images containing all elements")
println("generatefibres(j,savepp) generate j synthetic images containing fibres")
println("generatebackground(j,savepp), generaterings(j,savepp), generateclusters(j,savepp), generateclustersfibres(j,savepp) are also available.")
println("j defaults to 25 and savepp defaults to false. If savepp = true, the raw localisation coordinate files are saved too.")
println("j and savepp can be omitted to use the default values.")
println("----------")

function makeimages(files)
    data = Array{Float32,4}(undef,100,100,1,length(files))
    for k in axes(data,4)
        temp = CSV.read(files[k],DataFrame)
        if in("x [nm]", names(temp))
            rename!(temp, Dict("x [nm]" => "x[nm]", "y [nm]" => "y[nm]"))
        elseif in("x", names(temp))
            rename!(temp, Dict("x" => "x[nm]", "y" => "y[nm]"))
        end
        select!(temp,["x[nm]", "y[nm]"])
        temp = Matrix(temp)
        data[:,:,1,k] = bincounts(Hist2D((temp[:,1],temp[:,2]), binedges=(range(minimum(temp[:,1]), minimum(temp[:,1])+3000, length=101),range(minimum(temp[:,2]), minimum(temp[:,2])+3000, length=101))))
    end
    
    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley
    return data
end

function generaterings(j=25,savepp=false)
    data = Array{Float32,4}(undef,100,100,1,j)
    truth = Array{Float32,2}(undef,4,j)

    @showprogress Threads.@threads for i in axes(data,4)
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()

        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigmadist_c = 0

        initn = 0
        lendist = Bernoulli(0)
        densdist_f = Bernoulli(0)
        widthdist_f = Bernoulli(0)
        deldist = Bernoulli(0)
        sigdist = Bernoulli(0)
        bend = 0
        bezdeg = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp,"ring")
        truth[:,i] .= [numrings, radii, width_r, dens_r]
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:5
        for j in 1:5
            axs1[i,j] = Axis(fig1[i,j])
            hidedecorations!(axs1[i,j])
            hidespines!(axs1[i,j])
        end
    end

    for i in 1:minimum([size(data,4), 25])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
        text!(axs1[i],0,90,text = string("no. = ", maximum([0, round(Int,truth[1,i])])), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,80,text = string("radius = ", truth[2,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,70,text = string("width = ", 2*truth[3,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,60,text = string("dens. = ", truth[4,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
    end
    display(fig1)
end

function generateclusters(j=25,savepp=false)
    data = Array{Float32,4}(undef,100,100,1,j)
    truth = Array{Float32,2}(undef,3,j)

    @showprogress Threads.@threads for i in axes(data,4)
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()
        
        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        initn = 0
        lendist = Bernoulli(0)
        densdist_f = Bernoulli(0)
        widthdist_f = Bernoulli(0)
        deldist = Bernoulli(0)
        sigdist = Bernoulli(0)
        bend = 0
        bezdeg = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp,"clusters")
        truth[:,i] = [numnuc_c,dupenum_c,sigma_c]
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:5
        for j in 1:5
            axs1[i,j] = Axis(fig1[i,j])
            hidedecorations!(axs1[i,j])
            hidespines!(axs1[i,j])
        end
    end

    for i in 1:minimum([size(data,4), 25])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
        text!(axs1[i],0,90,text = string("no. = ", maximum([0, round(Int,truth[1,i])])), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,80,text = string("p.p.c. = ", truth[2,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,70,text = string("width = ", truth[3,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
    end

    display(fig1)
end

function generatebackground(j=25,savepp=false)
    data = Array{Float32,4}(undef,100,100,1,j)
    truth = Array{Float32,2}(undef,3,j)

    @showprogress Threads.@threads for i in axes(data,4)
        numnuc_b, dupedist_b, sigmadist_b = getbackground()

        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigmadist_c = Bernoulli(0)
        
        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        initn = 0
        lendist = Bernoulli(0)
        densdist_f = Bernoulli(0)
        widthdist_f = Bernoulli(0)
        deldist = Bernoulli(0)
        sigdist = Bernoulli(0)
        bend = 0
        bezdeg = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp,"background")
        truth[:,i] = [numnuc_b,dupenum_b,sigma_b]
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:5
        for j in 1:5
            axs1[i,j] = Axis(fig1[i,j])
            hidedecorations!(axs1[i,j])
            hidespines!(axs1[i,j])
        end
    end

    for i in 1:minimum([size(data,4), 25])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
        text!(axs1[i],0,90,text = string("no. = ", maximum([0, round(Int,truth[1,i])])), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,80,text = string("p.p.c. = ", truth[2,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,70,text = string("width = ", truth[3,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
    end

    display(fig1)
end

function generateall(j=25,savepp=false)
    data = Array{Float32,4}(undef,100,100,1,j)

    @showprogress Threads.@threads for i in axes(data,4)
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()

        data[:,:,1,i],numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp,"all")
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:5
        for j in 1:5
            axs1[i,j] = Axis(fig1[i,j])
            hidedecorations!(axs1[i,j])
            hidespines!(axs1[i,j])
        end
    end

    for i in 1:minimum([size(data,4), 25])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
    end
    display(fig1)
end

function generatefibres(j=25,savepp=false) 
    data = Array{Float32,4}(undef,100,100,1,j)
    truth = Array{Float32,2}(undef,5,j)

    @showprogress Threads.@threads for i in axes(data,4)
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        
        numnuc_b, dupedist_b, sigma_b = getbackground()

        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigma_c = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigma_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigma_b,savepp,"fibres")
        truth[:,i] = [numfibres, dens, width, curv, len]
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:25
        axs1[i] = Axis(fig1[fld1(i,5),mod1(i,5)])
        hidedecorations!(axs1[i])
        hidespines!(axs1[i])
    end

    for i in 1:minimum([size(data,4), 25])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
        text!(axs1[i],0,90,text = string("no. = ", maximum([0, round(Int,truth[1,i])])), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,80,text = string("dens. = ", truth[2,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,70,text = string("width = ", 2*truth[3,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,60,text = string("curv. = ", truth[4,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
        text!(axs1[i],0,50,text = string("length = ", truth[5,i]), color = :white, glowcolor = (:black, 1.0), glowwidth = 2.0)
    end
    

    display(fig1)
end

function generateclustersfibres(j=25,savepp=false)
    data = Array{Float32,4}(undef,100,100,1,j)
  
    @showprogress Threads.@threads for i in axes(data,4)
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()

        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        data[:,:,1,i], numfibres, dens, width, curv, len, delta, sigma, dupenum_c, sigma_c, numrings, radii, width_r, dens_r, dupenum_b, sigma_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp,"clustsfibres")
    end

    scaley = 40
    clamp!(data,0,scaley)
    data ./= scaley

    fig1 = Figure(size=(1000,1000))
    axs1 = Array{Any,2}(undef,5,5)
    for i in 1:25
        axs1[i] = Axis(fig1[fld1(i,5),mod1(i,5)])
        hidedecorations!(axs1[i])
        hidespines!(axs1[i])
        heatmap!(axs1[i],data[:,:,1,i],colorrange=[0,1])
    end
    display(fig1)
end
;