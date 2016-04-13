function [coreg,speg,keyg]=MICRFt(nodefile,texp,expgfile,netfile,outputfile,pi0,flagp)
%%
% Inputs:
% nodefile: node score file with two columns: Gene and score, tab separated.
%           When this is the only input, genes should be listed with GeneName. (Required)
% texp: time-series expressed data in this version, (dynamic network will
% included in future). (Required)
% expgfile: file with expressed genes. (Required)
% netfile:  network speific input by users with three clolumns gene1 gene2
%           and betweenness; one line with one edge. Or a specific string to select
%           one network used in MICRF paper, which should be one of the following
%           strings: 'STRING', 'iRefIndex', 'HPRD', 'CORR', 'CoEXP', 'CoPrePPI'. 
%           Default is 'CoPrePPI'. (Optional)
% outputfile:   output file name; if not given, it will only return a cell 
%               variable. (Optional)
% pi0:  the prior fraction of non-risk genes, default is 0.94. (Optional)
%
% Outputs:
% out:  a cell variable the same with output file
% w:    the optimal parameters used in MICRF model

%% initial inputs
netstrs={'STRING', 'iRefIndex', 'HPRD', 'CORR', 'CoEXP', 'CoPrePPI'};
netpath='data/network/';
netfiles={'STRINGnetmap.txt', 'iRefIndexm.txt', 'StringNew_HPRD_mnet.txt', 'brainspan_net_cor.txt', 'brainspan_net_top5.txt', 'ComCo_PrePPI.txt'};
bepath='data/Network_betweenness/';
benesss={'Betweenness_edge_STRING.txt', 'Betweenness_edge_iRef.txt', 'Betweenness_edge_HPRD.txt', 'Betweenness_edge_corr1.txt', 'Betweenness_edge_coexp5.txt', 'Betweenness_edge_Co_PrePPI.txt'};
Ws=zeros(2,6);
Ws(1,:)=[6.8557,6.1026,7.3602,7.4168,3.9194,3.2318];
Ws(2,:)=[10.653,9.4211,9.0906,8.0247,13.117,12.861];
flag=0; % 0 read MICRF files and 1 read users' input files

if nargin < 3
    error('At least two input parameters: nodefile and texp are needed for time-series MICRF\n');
end

if nargin == 3
    j=6;
    netfile=[netpath netfiles{j}];
    beness=[bepath benesss{j}];
    w0=Ws(:,j);
end

if nargin >= 4
    [~,j] = ismember(netfile,netstrs);
    if j > 0
        netfile=[netpath netfiles{j}];
        beness=[bepath benesss{j}];
        w0=Ws(:,j);
    elseif j==0
        flag = 1;
        w0=[mean(Ws(1,:));mean(Ws(2,:))];
    end
end

if nargin < 5
    outputfile='';
end

if nargin < 6
    pi0=0.94;
end

if nargin < 7
    flagp=0;
end
%% read related files
addpath(genpath(pwd)) % use functions in UGM package
% read the network edges and betweenness
if flag == 1
    [node1,node2,be] = textread(netfile,'%s%s%f','delimiter','\t');
    enode1=node1;enode2=node2;
elseif flag == 0
    [node1,node2,~] = textread(netfile,'%s%s%f','delimiter','\t');
    [enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
end
genes = union(node1,node2);
nStates = ones(1,length(genes)) * 2;
% read input node score
fcon =  fopen(nodefile,'r');
C = textscan(fcon,'%s%s','delimiter','\t');
fclose(fcon);
% read expressed genes
fcon =  fopen(expgfile,'r');
expG = textscan(fcon,'%s');
expG=expG{:};
fclose(fcon);

%% step 1: time-series node and edge features AND network adjacency matrix 
adj = step1_adjt(node1,node2,genes);
% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
[nodet,edget,t] = timefeas(texp,expG,genes,edgeStruct);

%% step 2 and 3: compute initial edge and node features; MICRF training, decoding and inferring
%  step 4: output files
outT=cell(t,1);
for i=1:t
    [Xnode,Xedge,nodeMap,edgeMap]=step2_featuret(genes,C,enode1,enode2,be,edgeStruct,pi0,nodet(:,i),edget(:,i));
    [Y,nps,w]=step3_MICRFt(Xnode,Xedge,nodeMap,edgeMap,edgeStruct,w0);
    if strcmp(outputfile,'') == 0
        outputfile1=[outputfile,int2str(i),'.txt'];
    else
        outputfile1=outputfile;
    end
    out=step4_outputt(genes,Y,nps,Xnode,outputfile1);
    outT{i}=out;
end

%% step 5: time-series risk genes output
%flagp=1; % whether plot the landscape figure
[coreg,speg,keyg]=step5_outputT(outT,texp,expG,t,outputfile,flagp);

end
