# Optimizer Test, Simulates CME fitted to Lotka-Volterra simulations (see SI Fig 1)

using Random, LinearAlgebra,DataFrames,Statistics,CSV,Plots

include("MaxEnt_Simplex.jl")
include("MaxEnt_Simplex_Fit.jl")

rng=MersenneTwister(1259);

AA=[0.6 1.2 4.0];

dt=.001;

# First value

it=1;

    A=[1 AA[it] AA[it]; AA[it] 1 AA[it]; AA[it] AA[it] 1];
    r=[1;1;1];
    numtim=3000000;

    x=ones(3,numtim+1);
    B=randn(rng,3,numtim);

    sigma=5.0;
    for t in 1:numtim
        x[:,t+1]=max.(x[:,t].+dt.*x[:,t].*r.*(ones(3,1).-A*x[:,t]).+sigma*dt.*B[:,t],0.0);
    end
    x=x[:,((x[1,:].>0).&(x[2,:].>0).&(x[3,:].>0))];
    x=x[:,30000:1000:2000000]
    x=[x[[1,2,3],:] x[[1,3,2],:] x[[2,1,3],:] x[[2,3,1],:] x[[3,2,1],:] x[[3,1,2],:]];

    x=x./sum(x,dims=1);

    lambda=-500;
    gamma=0;

    K1,h1,M1,M_t1,C1,C_t1=CME(x,lambda,gamma,1000000);


# Second value

it=2;

    A=[1 AA[it] AA[it]; AA[it] 1 AA[it]; AA[it] AA[it] 1];
    r=[1;1;1];
    numtim=3000000;

    x=ones(3,numtim+1);
    B=randn(rng,3,numtim);

    sigma=5.0;
    for t in 1:numtim
        x[:,t+1]=max.(x[:,t].+dt.*x[:,t].*r.*(ones(3,1).-A*x[:,t]).+sigma*dt.*B[:,t],0.0);
    end
    x=x[:,((x[1,:].>0).&(x[2,:].>0).&(x[3,:].>0))];
    x=x[:,30000:1000:2000000];
    x=[x[[1,2,3],:] x[[1,3,2],:] x[[2,1,3],:] x[[2,3,1],:] x[[3,2,1],:] x[[3,1,2],:]];


    x=x./sum(x,dims=1);

    lambda=-500;
    gamma=0;
    K2,h2,M2,M_t2,C2,C_t2=CME(x,lambda,gamma,1000000);


# Third value

  it=3;

    A=[1 AA[it] AA[it]; AA[it] 1 AA[it]; AA[it] AA[it] 1];
    r=[1;1;1];
    numtim=3000000;

    x=ones(3,numtim+1);
    B=randn(rng,3,numtim);

    sigma=5.0;
    for t in 1:numtim
        x[:,t+1]=max.(x[:,t].+dt.*x[:,t].*r.*(ones(3,1).-A*x[:,t]).+sigma*dt.*B[:,t],0.0);
    end
    x=x[:,((x[1,:].>0).&(x[2,:].>0).&(x[3,:].>0))];
    x=x[:,30000:1000:2000000];
    x=[x[[1,2,3],:] x[[1,3,2],:] x[[2,1,3],:] x[[2,3,1],:] x[[3,2,1],:] x[[3,1,2],:]];

    x=x./sum(x,dims=1);

    lambda=-500;
    gamma=0;
    K3,h3,M3,M_t3,C3,C_t3=CME(x,lambda,gamma,1000000);



# Plotting


plot(.31:.02:.365,.31:.02:.365,line=(:dot,5),color=:grey,label="",xtickfontsize=12,ytickfontsize=12,xguidefontsize=14,yguidefontsize=14,legendfontsize=14)
scatter!([M1 M2 M3],[M_t1 M_t2 M_t3],color=[1 2 3], xaxis= ("Sample Means", (.31,.365),.32:.02:.36),yaxis= ("Model Means", (.31,.365),.32:.02:.36),label=["\\alpha=0.6" "\\alpha=1.2" "\\alpha=4.0"],legend=:topleft,marker=7)


plot(-.05:.01:.10,-.05:.01:.10,line=(:dot,5),color=:grey,label="",xtickfontsize=12,ytickfontsize=12,xguidefontsize=14,yguidefontsize=14,legendfontsize=14)
scatter!([C1 C2 C3],[C_t1 C_t2 C_t3],color=[1 2 3], xaxis= ("Sample Covariances", (-.05,.11),-.05:.05:.10),yaxis= ("Model Covariances", (-.05,.11),-.05:.05:.10),label=["\\alpha=0.6 , K>0" "\\alpha=1.2 , K=0" "\\alpha=4.0 , K<0"],legend=:topleft,marker=7)
