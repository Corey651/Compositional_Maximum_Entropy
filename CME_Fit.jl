# Fits Compositional data to the CME model using L2 regularization gamma (gamma=0 in the associated manuscript)

using Statistics, Optim, SpecialFunctions, LinearAlgebra, LineSearches, Zygote

using IrrationalConstants:twoπ,halfπ,sqrtπ,sqrt2π,invπ,inv2π,invsqrt2,invsqrt2π,logtwo,logπ,log2π

import ChainRulesCore

import SpecialFunctions.logerfc
logerfc(x::Real) = _logerfc(float(x))

function _logerfc(x::Real)
        return log(erfcx(x)) - x^2
    else
        return log(erfc(x))
    end
end


ChainRulesCore.@scalar_rule(erf(x, y), (- (2 * exp(-x^2)) / sqrtπ, (2 * exp(-y^2)) / sqrtπ))
ChainRulesCore.@scalar_rule(logerfc(x), - (2 * exp(-x^2 - Ω)) / sqrtπ)
ChainRulesCore.@scalar_rule(logerf(x,y), (- (2 * exp(-x^2 - Ω)) / sqrtπ, (2 * exp(-y^2 - Ω)) / sqrtπ))


function CME_Fit(Data,gamma)

    # Inputs and Outputs

         #=

         Here the input Data is an NxD matrix
         N is the number of variables of interest, D is the number of samples
         The sum of the N variables (for each sample) needs to be normalized
         to 1


         The output Params is an N x N-1 matrix of maximum entropy parameters

         =#


    # Optimization

    #removing the redundant Nth variable

    sz=size(Data);

    Data=Data[1:sz[1]-1,:];

    M=mean(Data,dims=2);

    Chi=Data*Data';

    Chi=Chi/sz[2];

    Params=zeros(Float64,sz[1],sz[1]-1);

    Threads.@threads for i in 1:sz[1]-1

        L=copy(Data);
        L[i,:].=1.0;
        CC=copy([Chi[:,i];Chi[i,i]]);
        CC[i]=copy(M[i]);

        ic=zeros(Float64,1,sz[1]);
        V=Chi[i,i]-M[i]^2;
        ic[i]=M[i]/V-1/M[i];
        ic[end]=-1/(V);
        c=2 .-sum(L,dims=1);


        obj(Par)=LogPseudo(Par,CC,L,c,gamma,i);
        function g!(G,x)
            G.=obj'(x)
        end


        res=optimize(obj,g!, ic,LBFGS(; m=5, linesearch=BackTracking(order=3)));  #slow lin
        count=0;
        temp=Optim.minimizer(res)
        temp=temp[1]
        while isnan(temp)&&count<20
            count+=1;
            ic=ic+mean(ic)/5*(rand(Float64, 1,sz[1])-0.5*ones(Float64,1,sz[1]));
            res=optimize(obj,g!, ic,LBFGS(; m=5, linesearch=BackTracking(order=3)));
            temp=Optim.minimizer(res)
            temp=temp[1]
        end

        Params[:,i]=Optim.minimizer(res);
        Params[end,i]=2*Params[end,i];


    end


    P=zeros(Float64,sz[1],sz[1]);
    h=[diag(Params); 0];
	P[1:sz[1]-1,1:sz[1]-1]=Params[1:sz[1]-1,1:sz[1]-1];
	for i in 1:sz[1]-1
    	P[i,i]=Params[sz[1],i];
	end

	Params=zeros(Float64,sz[1],sz[1]);
	for i in 1:sz[1]
		for j in 1:sz[1]
			Params[i,j]=P[i,j]+P[j,i]-P[i,i]-P[j,j];
    	end
    end
    K=Params./2;
	h=h.-K[1:sz[1],sz[1]];

    return K,h

end

function LogPseudo(Param,Cons,Data,c,gamma,ii)

    PP=copy(Param[end]);
    PPP=copy([Param[1:1:ii-1]; Param[ii+1:1:end-1];0]);
    PPP=PPP.-2*copy(PP);               #to compensate for the ignored 1/2 term in the probability distribution

    if PP<0
        b=Param[1:1,1:end-1]*Data;
        c=2*sqrt(-PP)^2*c;
        c=b-c;
        c=c/(2*sqrt(-PP));
        b=b/(2*sqrt(-PP));
        c=logerf.(c,b);           #slow line
        D=length(c);
        return -1*dot(Param,Cons)-.5*log(-4*PP/pi)+sum(b.^2)/D+sum(c)/D+gamma*norm(PPP)^2
    else
        b=Param[1:1,1:end-1]*Data;
        c=2*sqrt(PP)^2*c;
        c=c+b;
        c=c/(2*sqrt(PP));
        b=b/(2*sqrt(PP));
        c=log.(erfi.(c)-erfi.(b));
        D=length(c);
        return -1*dot(Param,Cons)-.5*log(4*PP/pi)-sum(b.^2)/D+sum(c)/D+gamma*norm(PPP)^2
    end

end
