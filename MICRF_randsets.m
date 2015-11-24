function MICRF_randsets(kk,nrand)
addpath(genpath(pwd))
[netfiles,benesss,adjfiles]=getFiles();

wop0= zeros(2,6);
wop0(:,1) =[1;1];
wop0(:,2) =[6.9419;24.7418];
wop0(:,3) =[7.6585;25.9918];
wop0(:,4) =[5.4660;16.1872];
wop0(:,5) =[3.5793;27.2895];
wop0(:,6) =[1;1];

%% different for simulation sets
if nrand == 1    
    % rand set 1: j = 2:9, i=1:20 for each net
    nSim=160;
    inputpath='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset_1/';
    inputstr='rand1';
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult_5/MICRFresult_';
end

if nrand == 2
    % rand set 2: leave-one mutation, 121 gold risk genes
    nSim=121;
    inputpath='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/leaveone4/';
    inputstr='rand2_';
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/leaveone4result_5/MICRFresult_';
end

if nrand == 3 
    % rand set 3: recurrent mutation, 20 simulations for each net
    nSim=20;
    inputpath='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset4/';
    inputstr='part3_';
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult4_5/MICRFresult_';
end

if nrand == 4
    % rand set 4: 1911 controls samples, 1 for each net
    nSim=1;
    inputpath='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/control/';
    inputstr='control1911';
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/control/v5/MICRFresult_';
end

%%
netj=mod(kk,nSim);
if netj == 0 
	netj=nSim;
end

netflag=floor((kk-1)/nSim)+1;
if netflag == 1
    instr = 'CRF_input';
    idm=1;
else
    instr = 'hotnet_input'; 
    idm=0;
end

if nrand == 1
    % different simulated sets
    netj_i = floor((netj-1)/20)+2;
    netj_j = mod(netj,20);
    if netj_j == 0 
        netj_j=20;
    end
    nodef=[inputpath,instr,inputstr,int2str(netj_i),'_',int2str(netj_j),'.txt'];
    outputfile=[outputstr,int2str(netflag),'_',int2str(netj_i),'_',int2str(netj_j),'.txt'];
else
    nodef=[inputpath,instr,inputstr,int2str(netj),'.txt'];
    outputfile=[outputstr,int2str(netflag),'_',int2str(netj),'.txt'];
end

w0=wop0(:,netflag);
netfile=netfiles{netflag};
beness=benesss{netflag};
adjf=adjfiles{netflag};

%% step 1: network adjacency matrix 
%  adj = step1_adj(netfile);
load(adjf);

%% step 2 and 3: compute initial edge and node features; MICRF training, decoding and inferring
[genes,Xnode,Xedge,nodeMap,edgeMap,edgeStruct]=step2_feature(netfile,beness,nodef,adj);
[Y,nps]=step3_MICRF(Xnode,Xedge,nodeMap,edgeMap,edgeStruct,w0);

%% step 4: output files
step4_output(genes,Y,nps,outputfile,idm);

end
