# Compares CME (K) to the logit-normal model (K_ln)

using Statistics,LinearAlgebra
function Compare_to_logit(Data)
    K,h=CME_Fit(Data,0)
    sz=size(Data);
    CC=cov(Data');

    DD=log.(Data[1:sz[1]-1,:]./Data[sz[1]:sz[1],:]);
    C=cov(DD');
    Ko_ln=zeros(sz[1],sz[1])
    Ko_ln[1:sz[1]-1,1:sz[1]-1]=-pinv(C);
    K_ln=zeros(sz[1],sz[1])
    for i in 1:sz[1]
        for j in 1:sz[1]
            K_ln[i,j]=Ko_ln[i,j]+Ko_ln[j,i]-Ko_ln[i,i]-Ko_ln[j,j];
        end
    end
    return CC,K,K_ln
end
