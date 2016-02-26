function wop0=getOne()
addpath(genpath(pwd))
%% different for simulation sets
nSim=100;
outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/all/MICRFs_';
wop0= zeros(2,6);

w1V=ones(600,1);
w2V=ones(600,1);
fV =1000000 * ones(600,1);

for kk = 1:600
    netj=mod(kk,nSim);
    if netj == 0 
        netj=nSim;
    end
    netflag=floor((kk-1)/nSim)+1;
    outputfile=[outputstr,int2str(netflag),'_',int2str(netj),'.txt'];
    
    if exist(outputfile,'file') ==2
        fcon =  fopen(outputfile,'r');
        C = textscan(fcon,'%d\t%d\t%d\t%f\t%f\t%12.8f\t%d\n','delimiter','\t');
        fclose(fcon);
    
        w1V(kk) = C{4}(1);
        w2V(kk) = C{5}(1);
        fV(kk)= C{6}(1);
    end
end

for i = 1:6
    f100 = fV(((i-1)*100+1):(i*100));
    sub0 = find(f100==min(f100));
    sub = (i-1)*100 + sub0;
    wop0(1,i)=w1V(sub(1));
    wop0(2,i)=w2V(sub(1));
end

dlmwrite('wop_6.txt',wop0);

%% different for simulation sets
[netfiles,benesss,adjfiles,nodefiles]=getFiles();
outputfile = 'percent_TADA.txt';
noder = textread('data_hq/trueg.txt','%s');
subT=zeros(1000,1);

fcon =fopen(outputfile,'w');

for netflag = 1:6
if netflag==1
    nodef=nodefiles{1};
else
    nodef=nodefiles{2};   
end
netfile=netfiles{netflag};
beness=benesss{netflag};
adjf=adjfiles{netflag};
wop=wop0(:,netflag);
load(adjf);

[genes,Xnode,Xedge,nodeMap,edgeMap,edgeStruct]=step2_feature(netfile,beness,nodef,adj);
[nodePot,edgePot] = UGM_CRF_makePotentials(wop,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);

if netflag == 1
	fc =  fopen('/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/Fmap0121.txt','r');
	C = textscan(fc,'%s%s','delimiter','\t');
	id1=C{1};
	id2=C{2};
	fclose(fc);
    [Lia,j]=ismember(genes,id1);
    subs = find(Lia > 0);
    sub2 = find(~strcmp(id2(j(subs)),'NA'));
    genes(subs(sub2))=id2(j(subs(sub2)));
end


[~,j] = ismember(noder,genes);
subT(:,1)=j;

% decoding a CRF model
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
Y = int32(Y');
a=find(Y==1);
a1=length(intersect(a,subT(1:100,1)))/length(unique(subT(1:100,1)));
a2=length(intersect(a,subT(1:200,1)))/length(unique(subT(1:200,1)));
a3=length(intersect(a,subT(1:500,1)))/length(unique(subT(1:500,1)));

fprintf(fcon,'The number of risk genes: %d \n',sum(Y==1));
fprintf(fcon,'The 100 percent of risk genes: %f \n',a1);
fprintf(fcon,'The 200 percent of risk genes: %f \n',a2);
fprintf(fcon,'The 500 percent of risk genes: %f \n',a3);

end
fclose(fcon);

end

