# Fit and simulates the CME model of compositional data.

# lambda is a simulation parameter (with suggested values described in the included examples).
# gamma is an L2 regularization parameter
# t is the total length of the Monte Carlo simulation

using Statistics,LinearAlgebra,Distributions


function CME(Data,lambda,gamma,t)
    K,h=CME_Fit(Data,gamma)

    sz=length(h);
    Sim=zeros(sz,t);
    KK=-K./(2*lambda);
    hh=-h./(2*lambda);
    Sig=sqrt(-1/(2*lambda));
    s=(1/sz).*ones(sz,1);
    q=sum(s);
    for i in 1:sz
        M=hh[i]+dot(KK[i,:],s)+1-q+s[i];
        sam=rand(truncated(Normal(M,Sig),0,1))
        q=q-s[i,1]+sam;
        s[i,1]=sam;
    end
    for tt in 1:t
        for i in 1:sz
            M=hh[i]+dot(KK[i,:],s)+1-q+s[i];
            sam=rand(truncated(Normal(M,Sig),0,1))
            q=q-s[i,1]+sam;
            s[i,1]=sam;
        end
        for i in 1:sz
            M=hh[i]+dot(KK[i,:],s)+1-q+s[i];
            sam=rand(truncated(Normal(M,Sig),0,1))
            q=q-s[i,1]+sam;
            s[i,1]=sam;
        end

        Sim[:,tt]=s;
    end

    M=mean(Data,dims=2);
    M_t=mean(Sim,dims=2);
    C=cov(Data');
    C_t=cov(Sim');
    C=reshape(C,(sz*sz,1));
    C_t=reshape(C_t,(sz*sz,1));
    return K,h,M,M_t,C,C_t
end
