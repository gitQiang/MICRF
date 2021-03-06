function MICRF_randset4(kk)
addpath(genpath(pwd))

netfiles=cell(6,1);
netfiles{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/STRINGnetmap.txt';
netfiles{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/hotnet/iRefIndexm.txt'; 
netfiles{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/StringNew_HPRD_mnet.txt'; 
netfiles{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_cor.txt'; 
netfiles{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_top5.txt'; 
netfiles{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/ComCo_PrePPI.txt';

benesss=cell(6,1);
benesss{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_STRING.txt';
benesss{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_iRef.txt'; 
benesss{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_HPRD.txt'; 
benesss{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_corr1.txt'; 
benesss{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_coexp5.txt'; 
benesss{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_Co_PrePPI.txt';

wop0= zeros(2,6);
wop0(:,1) =[6;26];
wop0(:,2) =[6.9419;24.7418];
wop0(:,3) =[7.6585;25.9918];
wop0(:,4) =[5.4660;16.1872];
wop0(:,5) =[3.5793;27.2895];
wop0(:,6) =[6;26];

adjfiles=cell(6,1);
adjfiles{1}='adj_STRING.mat'; 
adjfiles{2}='adj_iRef.mat'; 
adjfiles{3}='adj_HPRD.mat'; 
adjfiles{4}='adj_braincor.mat'; 
adjfiles{5}='adj.mat'; 
adjfiles{6}='adj_Co_PrePPI.mat';  

netj=mod(kk,20);
if netj == 0 
	netj=20;
end
netflag=floor((kk-1)/20)+1;
if netflag == 1
instr = 'CRF_input';
idm=1;
else
instr = 'hotnet_input'; 
idm=0;
end

w0=wop0(:,netflag);
nodef=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset4/' instr 'part3_' int2str(netj) '.txt'];
%onefile=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset4/' instr 'part3_' int2str(netj) '.mat'];
outputfile=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult4_5/MICRFresult_' int2str(netflag) '_' int2str(netj) '.txt'];
netfile=netfiles{netflag};
beness=benesss{netflag};
adjf=adjfiles{netflag};

%% step 1: network adjacency matrix 
%  adj = step1_adj(netfile);
load(adjf);

%% step 2 and 3: compute initial edge and node features; MICRF training, decoding and inferring
[genes,Y,nps]=step3_MICRF(netfile,beness,nodef,adj,w0);

%% step 4: output files
step4_output(genes,Y,nps,outputfile,idm);

end
