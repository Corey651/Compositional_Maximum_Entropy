
Tab=readtable('Metabric_breast.xlsx');
fid=fopen('gene_list_288.txt');
Genes=textscan(fid,'%s');
fclose(fid);
Genes=string(Genes{1});

% First select the triple-negative patients

Tab.ERStatus=string(Tab.ERStatus);
Tab.PRStatus=string(Tab.PRStatus);
Tab.HER2Status=string(Tab.HER2Status);
ind=find((Tab.ERStatus=="Negative")&(Tab.PRStatus=="Negative")&(Tab.HER2Status=="Negative"));
TTab=Tab((Tab.ERStatus=="Negative")&(Tab.PRStatus=="Negative")&(Tab.HER2Status=="Negative"),:);

% Include the stationary (or "Invariant, Inv") measure information

Data_onco=readmatrix('Metabric_RNA_Inv_onco.csv');
Data_onco=Data_onco(:,ind);

Adj=readtable('adj288.csv');
Adj=Adj(:,2:289);

% Remove the CD79A immune cluster, all genes with .5 Correlation or higher

CCC=corrcoef(Data_onco');
indd=find(CCC(:,41)>=.50);
ind=setdiff(1:288,indd);
Genes=Genes(ind);
Data_onco=Data_onco(ind,:);
Adj=Adj(ind,ind);


CC=cov(Data_onco');
DD=diag(CC);
[a,b]=sort(DD);
scatter(1:length(a),a)   % find the kink at 254
ind=b(254:end);
Genes=Genes(ind);

Data_onco=Data_onco(ind,:);
Adj=Adj(ind,ind);
Adj=table2array(Adj);
Adj=Adj-diag(diag(Adj));

CCCC=corrcoef(Data_onco');

% Reorder the genes to reveal covariance structures

[aa,bb]=Clusterfunc(CCCC,10);
Genes=Genes(aa);
Data_onco=Data_onco(aa,:);
Adj=Adj(aa,aa);

%writematrix(Data_onco,'Metabric_Input.csv')
%writematrix(Genes,'Genes_Filtered.csv')
%writematrix(Adj,'Adj_Filtered.csv')



