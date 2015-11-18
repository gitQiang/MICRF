function MICRF_randset4(kk)

addpath(genpath(pwd))
netfiles=cell(6,1);
netfiles{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/STRINGnetmap.txt';
netfiles{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/hotnet/iRefIndexm.txt'; netfiles{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/StringNew_HPRD_mnet.txt'; netfiles{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_cor.txt'; netfiles{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_top5.txt'; netfiles{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/ComCo_PrePPI.txt';
benesss=cell(6,1);
benesss{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_STRING.txt';
benesss{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_iRef.txt'; benesss{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_HPRD.txt'; benesss{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_corr1.txt'; benesss{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_coexp5.txt'; benesss{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_Co_PrePPI.txt';
wop0= zeros(2,6);
wop0(:,1) =[];
wop0(:,2) =[6.9419;24.7418];
wop0(:,3) =[];
wop0(:,4) =[];
wop0(:,5) =[3.5793;27.2895];
wop0(:,6) =[1;1];

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
w0=zeros(2,1);
w0=wop0(:,netflag);
nodef=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset4/' instr 'part3_' int2str(netj) '.txt'];
onefile=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset4/' instr 'part3_' int2str(netj) '.mat'];
outputfile=['/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult4_5/MICRFresult_' int2str(netflag) '_' int2str(netj) '.txt'];
netfile=netfiles{netflag};
beness=benesss{netflag};
step4_output(netfile,beness,nodef,netflag,onefile,outputfile,w0,idm)

end
