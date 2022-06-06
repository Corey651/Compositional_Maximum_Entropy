using DelimitedFiles, Distributed, Plots, Statistics, Colors

Data=readdlm("Metabric_Input.csv",',',Float64);
Genes=readdlm("Genes_Filtered.csv",',',String)
Adj=readdlm("Adj_Filtered.csv",',',Float64);


include("CME_Fit.jl")
include("Compare_to_logit.jl")
addprocs(8)
Data=Data./sum(Data,dims=1)

C,K,K_ln=Compare_to_logit(Data);

rmprocs(8)
sz=size(K);
KK=K.-sum(sum(K))/(sz[1]*(sz[1]-1));
KK=KK-Diagonal(KK);

C=C-Diagonal(C);

KK_ln=K_ln.-sum(sum(K_ln))/(sz[1]*(sz[1]-1));
KK_ln=KK_ln-Diagonal(KK_ln);

pyplot()

heatmap(KK,xtickfontrotation=90,xtickfontsize=12,xticks=(1:sz[1],vec(Genes)),ytickfontsize=12,yticks=(1:sz[1],vec(Genes)),c=reverse(colormap("RdBu",50,logscale=true)),clim=(-60000,60000),colorbar_title="Interaction",colorbar_ticks=[],colorbar_titlefontsize=14)
#savefig("Fig4a.pdf")
heatmap(C,clim=(-1.0e-5,1.0e-5),xtickfontsize=12,xtickfontrotation=90,xticks=(1:sz[1],vec(Genes)),ytickfontsize=12,yticks=(1:sz[1],vec(Genes)),colorbar_title="Covariance",colorbar_ticks=[],c=reverse(colormap("RdBu",50,logscale=true)),colorbar_titlefontsize=14)
#savefig("Fig4c.pdf")

heatmap(KK_ln,clim=(-500,500),xtickfontsize=12,xtickfontrotation=90,xticks=(1:sz[1],vec(Genes)),ytickfontsize=12,yticks=(1:sz[1],vec(Genes)),colorbar_title="Interaction",colorbar_ticks=[],c=reverse(colormap("RdBu",50,logscale=true)),colorbar_titlefontsize=14)
#savefig("Fig4b.pdf")
