function MICRFs(kk)
addpath(genpath(pwd))
[netfiles,benesss,adjfiles,nodefiles]=getFiles();

%% different for simulation sets
nSim=100;
outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/all/MICRFs_';

netj=mod(kk,nSim);
if netj == 0 
	netj=nSim;
end
netflag=floor((kk-1)/nSim)+1;
outputfile=[outputstr,int2str(netflag),'_',int2str(netj),'.txt'];

if netflag==1
    nodef=nodefiles{1};
else
    nodef=nodefiles{2};   
end
netfile=netfiles{netflag};
beness=benesss{netflag};
adjf=adjfiles{netflag};

load(adjf);
[~,Xnode,Xedge,nodeMap,edgeMap,edgeStruct]=step2_feature(netfile,beness,nodef,adj);
oneLearn(netj,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,outputfile);

end
