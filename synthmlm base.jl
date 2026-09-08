import Base.rand

# DO NOT RUN THIS FILE

# This file contains functions needed for SynthMLM and provides no functionality alone.

function getfibreset()
    initn = rand(1:50)

    avglen = rand(Uniform(1000,6000))
    lenrange = rand(Uniform(0.1,0.5))
    lendist = truncated(Normal(avglen,lenrange*avglen),lower = 500, upper = 10000)

    avgdens = rand(Uniform(0.8,2))
    densrange = rand(Uniform(0.1,0.5))
    densdist = truncated(Normal(avgdens,densrange*avgdens),lower = 0.25, upper = 6)

    fibwidth = rand(Uniform(10,100))
    widthrange = rand(Uniform(0.1,0.5))
    widthdist = truncated(Normal(fibwidth,widthrange*fibwidth),lower = 5, upper=300)

    avgdelta = rand(Uniform(0.025*avglen,0.25*avglen))
    deltarange = rand(Uniform(0.1,0.5))
    deldist = truncated(Normal(avgdelta,deltarange*avgdelta),lower=0.01*avglen,upper=0.5*avglen)

    avgsig = rand(Uniform(0.25*avgdelta,avgdelta))
    sigrange = rand(Uniform(0.1,0.5))
    sigdist = truncated(Normal(avgsig,sigrange*avgsig),lower=0.5*avgsig,upper=3*avgdelta)

    bend = rand(Uniform(0.2,1))
    bezdeg = rand(2:5)

    return initn, lendist, densdist, widthdist, deldist, sigdist, bend, bezdeg
end

function getringset()
    numrings = rand(1:400)

    avgrad = rand(Uniform(55,120))
    radrange = rand(Uniform(0.1,0.5))
    raddist = truncated(Normal(avgrad,radrange*avgrad),lower=25,upper=300)

    widupper = minimum([20, 0.7*avgrad/3])
    avgwidth = rand(Uniform(5,widupper))
    widrange = rand(Uniform(0.1,0.5))
    widthdist = truncated(Normal(avgwidth,widrange*avgwidth),lower=3,upper=3*widupper)

    avgdens = rand(Uniform(0.2,1.5))
    densrange = rand(Uniform(0.1,0.5))
    densdist = truncated(Normal(avgdens,densrange*avgdens),lower=0.1,upper=4.5)

    return numrings, raddist, widthdist, densdist
end

function getclustset()
    numnuc_c = rand(1:800)

    avgsig = rand(Uniform(5,100))
    sigrange = rand(Uniform(0.1,0.5))
    sigmadist_c = truncated(Normal(avgsig,sigrange*avgsig),lower=5,upper=300)

    maxdupe = 8*avgsig
    avgdupe =  rand(11:minimum([maxdupe, 400]))
    duperange = rand(Uniform(0.1,0.5))
    dupedist_c = truncated(Normal(avgdupe,duperange*avgdupe),lower=11,upper=minimum([maxdupe, 1000]))

    return numnuc_c, dupedist_c, sigmadist_c
end

function getbackground()
    numnuc_b = rand(1:500)

    avgdupe =  rand(1:10)
    duperange = rand(Uniform(0.1,0.5))
    dupedist_b = truncated(Normal(avgdupe,duperange*avgdupe),lower=1,upper=10)

    avgsig = rand(Uniform(5,60))
    sigrange = rand(Uniform(0.1,0.5))
    sigmadist_b = truncated(Normal(avgsig,sigrange*avgsig),lower=5,upper=180)

    return numnuc_b, dupedist_b, sigmadist_b
end

function bez(t,P)
    bez = @SVector[.0,.0]
    n = size(P,1)-1

    for i in axes(P,1)
        bez += bezcoef(i-1,n,t) .* P[i,:]
    end

    return bez
end

function phi(x)
    return 1 / (sqrt(2 * pi)) * exp(-x^2/2)
end

function CDF(x,delta,sig,N,fuzz)
    CDF = .0
    for i in -N:2*N
        CDF += (1 + erf((x - i*delta-fuzz[i+N+1])/(sig * sqrt(2)))) / (2+(3*N+1))
    end

    return CDF
end

function getbenddeg(curv,curvs,bends,bezdegs)
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

    if bezdegs > 3
        bezdegs -= 1
    end

    return bends, bezdegs
end

function getCU(len,delta,sig,np)
    if delta == 0
        delta = 0.1
    end
    
    N = floor(Int,len/delta)
    fuzz = 0.2*delta*rand(truncated(Normal(), -delta/2, delta/2),(3*N+1)) .+ rand(Uniform(-delta/2,delta/2))
    cee = CDF.(range(0,len,2*np),delta,sig,N,Ref(fuzz))
    if cee[1] == cee[end]
        U = fill(cee[1],np)
    else
        U = rand(Uniform(cee[1],cee[end]),np)
    end
    sort!(U)

    return cee, U
end

function getts(cee,U,len,np,P)
    ts = Vector{Float64}(undef,np)
    lenmod = len/length(cee)
    lenints = 1000
    ints = Vector{Float64}(undef,lenints)
    ints[1] = 0
    for i in 2:lenints
        ints[i] = ints[i-1] + arclength((i-1)/lenints,i/lenints,P)
    end
    for i in eachindex(ts)
        ts[i] = (searchsortedfirst(ints,(searchsortedfirst(cee, U[i])-1)*lenmod)-1) / lenints
    end

    return ts
end

function gettlim(t,bezt)
    inbox = bezt[:,1] .>=0 .&& bezt[:,1] .<= 3000 .&& bezt[:,2] .>= 0 .&& bezt[:,2] .<= 3000
    lims = Vector{Union{Int, Nothing}}(undef,2)
    lims[1] = findfirst(inbox)
    if isnothing(lims[1])
        return bezt[1,:], t[1]
    end

    lims[2] = findfirst(.!inbox[lims[1]:end])
    if isnothing(lims[2])
        lims[2] = length(t)
    else
        lims[2] += lims[1] - 1
    end

    sort!(lims)
    
    return bezt[lims[1]:lims[2],:], t[lims[1]:lims[2]]
end

function bezandchop(cee,U,len,np,P)
    protots = getts(cee,U,len,np,P)
    protobezt = Array{Float64,2}(undef,length(protots),2)

    for i in eachindex(protots)
        protobezt[i,:] .= bez(protots[i],P)
    end

    bezt, ts = gettlim(protots,protobezt)
    return bezt, ts
end

function bezandchopOG(ss,len,np,P)
    protots = Vector{Float64}(undef,np)
    protots[1] = ss[1] / norm(bezdiv(0,P))
    for i in 2:np
        protots[i] = protots[i-1] + (ss[i] − ss[i-1]) / norm(bezdiv(protots[i-1],P)) # rough and ready.
    end

    protobezt = Array{Float64,2}(undef,length(protots),2)

    for i in eachindex(protots)
        protobezt[i,:] .= bez(protots[i],P)
    end

    bezt, ts = gettlim(protots,protobezt)
    return bezt, ts
end


function bezpointsandfacts(P,dens,width,delta,sig)
    len = arclength(0,1,P)
    np = round(Int,dens*len)
    if np < 1
        np = 1
    end
    cee, U = getCU(len,delta,sig,np)
    
    bezt, ts = bezandchop(cee,U,len,np,P)

    len2 = arclength(ts[1],ts[end],P)
    np2 = length(ts)
    pp = Array{Float64,2}(undef,np2,2)
    widths = width*randn(round(Int,np2))

    for i in axes(pp,1)
        pp[i,:] = bezt[i,:] .+ widths[i] * normvec(ts[i],P)
    end

    curv = mean(abs.(bezcurv.(ts[1]:0.001:ts[end],Ref(P)))) # before was bias
    # curv = mean(abs.(bezcurv.(ts,Ref(P))))
    return pp, curv, len2
end

function bezpointsandfactsOG(P,dens,width,delta,sig)
    len = arclength(0,1,P)
    np = round(Int,dens*len)
    
    ss = sort!(rand(Uniform(0,len),np))
    
    bezt, ts = bezandchopOG(ss,len,np,P)

    len = abs(arclength(ts[1],ts[end],P))

    np2 = length(ts)
    pp = Array{Float64,2}(undef,np2,2)
    widths = width*randn(round(Int,np2))

    for i in axes(pp,1)
        pp[i,:] = bezt[i,:] .+ widths[i] * normvec(ts[i],P)
    end

    curv = mean(abs.(bezcurv.(ts,Ref(P))))

    return pp, curv, len
end

function getP(bezdeg,projlen,bend,phi,shift)
    P = Array{Float64,2}(undef,bezdeg+1,2)

    for i in 2:bezdeg
        P[i,:] .= [cos(phi) -sin(phi); sin(phi) cos(phi)] * rand(Uniform(-bend*projlen,bend*projlen),2) .+ shift
    end
    P[1,:] .= [cos(phi) -sin(phi); sin(phi) cos(phi)] * [-projlen/2, 0] .+ shift
    P[end,:] .= [cos(phi) -sin(phi); sin(phi) cos(phi)] * [projlen/2, 0] .+ shift

    return SMatrix{bezdeg+1,2}(P)
end

function getbez(projlen,dens,width,bend,bezdeg,phi,shift,delta,sig)
    P = getP(bezdeg,projlen,bend,phi,shift)
    pp, meancurv, length = bezpointsandfacts(P,dens,width,delta,sig)
    
    return pp, meancurv, length
end

function getbezOG(projlen,dens,width,bend,bezdeg,phi,shift,delta,sig)
    P = getP(bezdeg,projlen,bend,phi,shift)

    pp, meancurv, length = bezpointsandfactsOG(P,dens,width,delta,sig)
    
    return pp, meancurv, length
end

function bezcoef(i,n,t)
    return binomial(n,i) * t^i * (1-t)^(n-i)
end

function bezdiv(t,P)
    Bprime = @SVector [.0,.0]
    n = size(P,1)-1
    for i in 1:n
        Bprime += n * bezcoef(i-1,n-1,t) * (P[i+1,:] - P[i,:])
    end

    return Bprime
end

function bezdivdiv(t,P)
    Bprimeprime = @SVector [.0,.0]
    n = size(P,1)-1
    for i in 1:(n-1)
        Bprimeprime += n * (n - 1) * bezcoef(i-1,n-2,t) * (P[i+2,:] -2*P[i+1,:] - P[i,:]) # should I have divided the curvature by 3000
    end

    return Bprimeprime
end

function bezcurv(t,P)
    return norm(bezdivdiv(t,P)) /(10^4)
end

function normvec(t,P)
    n = bezdivdiv(t,P)./norm(bezdiv(t,P))^2 .- dot(bezdivdiv(t,P)./norm(bezdiv(t,P))^2,bezdiv(t,P)./norm(bezdiv(t,P))) .* bezdiv(t,P)./norm(bezdiv(t,P))
    return n ./ norm(n)
end

function arclength(t1,t2,P)
    f(u,p) = norm(bezdiv(u,P)) # p is a parameter vector, which there is none, and P is the control points, of which there are some.
    domain = (t1, t2)
    prob = IntegralProblem(f, domain)
    sol = Integrals.solve(prob, HCubatureJL())
    
    return sol.u
end

function makefibre(lendist,densdist,widthdist,deldist,sigdist,bend,bezdeg)
    phi = rand(Uniform(0,2*pi))
    shift = rand(Uniform(0,3000),2)
    projlen = rand(lendist)
    dens = rand(densdist) 
    width = rand(widthdist)
    delta = rand(deldist)
    sigma = rand(sigdist)
    
    if projlen != 0
        pp, curv, length = getbez(projlen,dens,width,bend,bezdeg,phi,shift,delta,sigma)
    else
        pp = []
        dens = 0
        width = 0
        curv = 0
        length = 0
        delta = 0
        sigma = 0
    end

    return pp, dens, width, curv, length, delta, sigma
end

function addclusters(numnuc,dupedist,sigmadist,pp = [])
    initpos = rand(Uniform(0,3000), numnuc,2)
    dupenum = round.(Int,rand(dupedist, numnuc))
    sigma = rand(sigmadist, numnuc)
    cluzzies = Array{Float64,2}(undef, sum(dupenum),2)
    for i in eachindex(dupenum)
        if any(dupenum .<0)
            println(dupenum)
        end
        prevnum = sum(dupenum[1:i-1])
        for j in prevnum+1:prevnum+dupenum[i]
            cluzzies[j,:] = initpos[i,:] .+ sigma[i] .* randn(2)
        end
    end

    if isempty(pp)
        pp = cluzzies
    else
        pp = [pp; cluzzies]
    end
    
    if numnuc > 1
        dupem = mean(dupenum)
        dupestd = std(dupenum,mean=dupem)

        sigm = mean(sigma)
        sigstd = std(sigma,mean=sigm)

        return pp, dupem, dupestd, sigm, sigstd
    elseif numnuc == 1
        dupem = mean(dupenum)

        sigm = mean(sigma)

        return pp, dupem, 0, sigm, 0
    else
        return pp, 0, 0, 0, 0
    end
end

function addsnow(numsnow,pp = [])
    initpos = rand(Uniform(0,3000), numnuc,2)
    
    sigma = rand(sigmadist, numnuc)
    cluzzies = Array{Float64,2}(undef, sum(dupenum),2)
    for i in eachindex(dupenum)
        if any(dupenum .<0)
            println(dupenum)
        end
        prevnum = sum(dupenum[1:i-1])
        for j in prevnum+1:prevnum+dupenum[i]
            cluzzies[j,:] = initpos[i,:] .+ sigma[i] .* randn(2)
        end
    end

    if isempty(pp)
        pp = cluzzies
    else
        pp = [pp; cluzzies]
    end
    
    if numnuc > 1
        dupem = mean(dupenum)
        dupestd = std(dupenum,mean=dupem)

        sigm = mean(sigma)
        sigstd = std(sigma,mean=sigm)

        return pp, dupem, dupestd, sigm, sigstd
    elseif numnuc == 1
        dupem = mean(dupenum)

        sigm = mean(sigma)

        return pp, dupem, 0, sigm, 0
    else
        return pp, 0, 0, 0, 0
    end
end

function addfibres(initn,lendist,densdist,widthdist,deldist,sigdist,bend,bezdeg,pp=[])
    dens = Vector{Float64}(undef,initn)
    width = Vector{Float64}(undef,initn)
    curv = Vector{Float64}(undef,initn)
    len = Vector{Float64}(undef,initn)
    deltas = Vector{Float64}(undef,initn)
    sigs = Vector{Float64}(undef,initn)

    pp0, dens[1], width[1], curv[1], len[1], deltas[1], sigs[1] = makefibre(lendist,densdist,widthdist,deldist,sigdist,bend,bezdeg)

    if isempty(pp)
        pp = pp0
    else
        pp = [pp; pp0]
    end

    for i in 2:initn
        pp_new, dens[i], width[i], curv[i], len[i], deltas[i], sigs[i] = makefibre(lendist,densdist,widthdist,deldist,sigdist,bend,bezdeg)
        pp = [pp; pp_new]
    end

    visfib =  len .> 0 .&& len .> 2 .* width
    numfibres = count(visfib)
    
    if numfibres > 1
        densm = mean(dens[visfib])
        denstd = std(dens[visfib],mean=densm)
        
        widm = mean(width[visfib])
        widstd = std(width[visfib],mean=widm)

        curvm = mean(curv[visfib])
        curvstd = std(curv[visfib],mean=curvm)

        lenm = mean(len[visfib])
        lenstd = std(len[visfib],mean=lenm)

        delm = mean(deltas[visfib])
        delstd = std(deltas[visfib],mean=delm)
        
        sigm = mean(sigs[visfib])
        sigstd = std(sigs[visfib],mean=sigm)

        return pp, numfibres, densm, denstd, widm, widstd, curvm, curvstd, lenm, lenstd, delm, delstd, sigm, sigstd
    elseif numfibres == 1
        densm = mean(dens[visfib])

        widm = mean(width[visfib])

        curvm = mean(curv[visfib])

        lenm = mean(len[visfib])

        delm = mean(deltas[visfib])
        
        sigm = mean(sigs[visfib])
        
        return pp, numfibres, densm, 0, widm, 0, curvm, 0, lenm, 0, delm, 0, sigm, 0
    else
        return pp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    end
end

function addrings(numrings,raddist,widthdist,densdist,pp=[])
    radii = rand(raddist,numrings)
    width = rand(widthdist,numrings)
    dens = rand(densdist,numrings) # density per unit length along circle mindens,maxdens

    centres = rand(Uniform(0,3000),numrings,2)
    shifty!(radii,centres)

    np = round.(Int, 2 * pi .* dens .* radii)
    xs = Vector{Float64}(undef,sum(np))
    ys = Vector{Float64}(undef,sum(np))

    for i in eachindex(np)
        prevnum = sum(np[1:i-1])
        mu = rand(Uniform(0,2*pi))
        K = 2*rand(Beta(2,8))
        for j in prevnum+1:prevnum+np[i]
            phi = rand(VonMises(mu,K)) # rand(Uniform(0,2*pi)) 
            wig = width[i]*randn()
            xs[j] = (wig + radii[i]) * cos(phi) + centres[i,1]
            ys[j] = (wig + radii[i]) * sin(phi) + centres[i,2]
        end
    end  
    
    pp0 = [xs ys]

    if isempty(pp)
        pp = pp0
    else
        pp = [pp; pp0]
    end

    hiddenring = centres[:,1] .<0 .|| centres[:,1] .> 3000 .|| centres[:,2] .< 0 .|| centres[:,2] .> 3000

    visrings = count(.!hiddenring)

    if visrings > 1
        radm = mean(radii[.!hiddenring])
        radstd = std(radii[.!hiddenring],mean=radm)

        widm = mean(width[.!hiddenring])
        widstd = std(width[.!hiddenring],mean=widm)

        densm = mean(dens[.!hiddenring])
        densstd = std(dens[.!hiddenring],mean=densm)

        return pp, visrings, radm, radstd, widm, widstd, densm, densstd
    elseif visrings == 1
        radm = mean(radii[.!hiddenring])

        widm = mean(width[.!hiddenring])

        densm = mean(dens[.!hiddenring])
        
        return pp, visrings, radm, 0, widm, 0, densm, 0
    else
        return pp, 0, 0, 0, 0, 0, 0, 0
    end
end

function shifty!(radii,centres)
    intercept = getintercepts(centres,radii)
    # needs to be a matrix
    while any(intercept) == true
        inds = findall(intercept)
        for a in inds
            i = a[1]
            j = a[2]
            diff = centres[j,:] - centres[i,:]
            diffdotdiff = diff[1]^2 + diff[2]^2
            sumrad = radii[i] + radii[j]
            lambda = (sqrt((sumrad^2+diffdotdiff)/diffdotdiff) - 1)/2
            centres[i,:] .-= lambda * diff
            centres[j,:] .+= lambda * diff
        end
        intercept = getintercepts(centres,radii)
    end

    return nothing
end

function getintercepts(centres,radii)
    intercept = fill(false, length(radii),length(radii))
    for i in axes(intercept,1)
        for j in i+1:length(radii)
            intercept[i,j] = norm(centres[i,:] - centres[j,:]) < radii[i]+ radii[j]
        end
    end
    return intercept
end

function makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b,savepp=false,name=[])
    pp, dupenum_b, dupestd_b, sigma_b, sigstd_b = addclusters(numnuc_b,dupedist_b,sigmadist_b)
    # snownuc = rand(0:500) # what about this!!! (500,2000) originally
    # pp, _, _, _, _ = addclusters(snownuc,1:1,30,pp)

    if initn > 0
        pp, numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd = addfibres(initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg,pp)
    else
        numfibres = 0
        dens = 0
        densstd = 0
        width = 0
        widthstd = 0
        curv = 0
        curvstd = 0
        len = 0
        lenstd = 0
        delta = 0
        deltastd = 0
        sigma = 0
        sigmastd = 0
    end

    if initnumrings > 0
        pp, numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r = addrings(initnumrings,raddist,widthdist_r,densdist_r,pp)
    else
        numrings = 0
        radii = 0
        radstd = 0
        width_r = 0
        widthstd_r = 0
        dens_r = 0
        densstd_r = 0
    end

    if numnuc_c > 0
        pp, dupenum_c, dupestd_c, sigma_c, sigstd_c = addclusters(numnuc_c,dupedist_c,sigmadist_c,pp)
    else
        dupenum_c = 0
        dupestd_c = 0
        sigma_c = 0
        sigstd_c = 0
    end
    
    binpp = bincounts(Hist2D((pp[:,1],pp[:,2]), binedges=((range(0,3000,length=101)),range(0,3000,length=101))))

    if savepp
        pp = pp[ (0 .<= pp[:,1] .< 3000) .& (0 .<= pp[:,2] .< 3000) ,:]
        
        CSV.write(string(name,"_synth_pp_",i,".csv"), DataFrame(pp,["x[nm]", "y[nm]"]))
        
        pp[:,1] .+= mod(i-1,10)*3000
        pp[:,2] .+= fld(i-1,10)*3000
        if i == 1
            appendy = false
        else
            appendy = true
        end
    end

    return binpp, numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b
end

function megaset(saveloc)
    nimages = 50000
    data = Array{Float32,4}(undef,100,100,1,nimages)
    truth = Array{Any,2}(undef,31,nimages)

    clustinds = 1:7500
    backinds = 7001:10000
    fibinds = 10001:20000
    fibandclustinds = 20001:30000
    ringinds = 30001:35000
    ringandclustinds = 35001:40000
    ringsandfibreinds = 40001:45000
    allinds = 45001:50000

    println("just background")
    @showprogress Threads.@threads for i in backinds
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

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,0, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b,"background"]
    end

    println("all")
    @showprogress Threads.@threads for i in allinds
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [numfibres,dens,densstd,width,widthstd,curv,curvstd,len,lenstd,delta,deltastd,sigma,sigmastd, numnuc_c,dupenum_c,dupestd_c,sigma_c,sigstd_c, numrings,radii,radstd,width_r,widthstd_r,dens_r,densstd_r, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "all"]
    end

    println("fibres and rings")
    @showprogress Threads.@threads for i in ringsandfibreinds
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        
        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigmadist_c = Bernoulli(0)

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [numfibres,dens,densstd,width,widthstd,curv,curvstd,len,lenstd,delta,deltastd,sigma,sigmastd, 0,0,0,0,0, numrings,radii,radstd,width_r,widthstd_r,dens_r,densstd_r, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "fibresandrings"]
    end

    println("fibres and clusters")
    @showprogress Threads.@threads for i in fibandclustinds
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()

        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [numfibres,dens,densstd,width,widthstd,curv,curvstd,len,lenstd,delta,deltastd,sigma,sigmastd, numnuc_c,dupenum_c,dupestd_c,sigma_c,sigstd_c, 0,0,0,0,0,0,0, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "fibresandclusts"] 
    end

    println("rings and clusters")
    @showprogress Threads.@threads for i in ringandclustinds
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        numnuc_c, dupedist_c, sigmadist_c = getclustset()
        
        initn = 0
        lendist = Bernoulli(0)
        densdist_f = Bernoulli(0)
        widthdist_f = Bernoulli(0)
        deldist = Bernoulli(0)
        sigdist = Bernoulli(0)
        bend = 0
        bezdeg = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [0,0,0,0,0,0,0,0,0,0,0,0,0, numnuc_c,dupenum_c,dupestd_c,sigma_c,sigstd_c, numrings,radii,radstd,width_r,widthstd_r,dens_r,densstd_r, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "ringsandclusts"]
    end

    println("clusters")
    @showprogress Threads.@threads for i in clustinds
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

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [0,0,0,0,0,0,0,0,0,0,0,0,0, numnuc_c,dupenum_c,dupestd_c,sigma_c,sigstd_c, 0,0,0,0,0,0,0, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "clusts"] 
    end

    println("fibres")
    @showprogress Threads.@threads for i in fibinds # 200k just fibres
        initn, lendist, densdist_f, widthdist_f, deldist, sigdist, bend, bezdeg = getfibreset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()
        
        initnumrings = 0
        raddist = Bernoulli(0)
        widthdist_r = Bernoulli(0)
        densdist_r = Bernoulli(0)

        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigmadist_c = Bernoulli(0)

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [numfibres,dens,densstd,width,widthstd,curv,curvstd,len,lenstd,delta,deltastd,sigma,sigmastd, 0,0,0,0,0, 0,0,0,0,0,0,0, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "fibres"] 
    end

    println("rings")
    @showprogress Threads.@threads for i in ringinds
        initnumrings, raddist, widthdist_r, densdist_r = getringset()
        numnuc_b, dupedist_b, sigmadist_b = getbackground()

        numnuc_c = 0
        dupedist_c = Bernoulli(0)
        sigmadist_c = Bernoulli(0)

        initn = 0
        lendist = Bernoulli(0)
        densdist_f = Bernoulli(0)
        widthdist_f = Bernoulli(0)
        deldist = Bernoulli(0)
        sigdist = Bernoulli(0)
        bend = 0
        bezdeg = 0

        data[:,:,1,i], numfibres, dens, densstd, width, widthstd, curv, curvstd, len, lenstd, delta, deltastd, sigma, sigmastd, dupenum_c, dupestd_c, sigma_c, sigstd_c,  numrings, radii, radstd, width_r, widthstd_r, dens_r, densstd_r, dupenum_b, dupestd_b, sigma_b, sigstd_b = makesynth(i,initn,lendist,densdist_f,widthdist_f,deldist,sigdist,bend,bezdeg, numnuc_c,dupedist_c,sigmadist_c, initnumrings,raddist,widthdist_r,densdist_r, numnuc_b,dupedist_b,sigmadist_b)
        truth[:,i] = [0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0, numrings,radii,radstd,width_r,widthstd_r,dens_r,densstd_r, numnuc_b,dupenum_b,dupestd_b,sigma_b,sigstd_b, "rings"]
    end

    @save saveloc data truth
end