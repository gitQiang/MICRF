fM = -1 * ones(10,10);
w1M = -1 * ones(10,10);
w2M = -1 * ones(10,10);

YM = -1 * ones(10,10);

for kk=1:100

    i = floor((kk-1)/10) + 1;
    j = mod(kk,10);
    if(j==0) j = 10; end
    
    filename= [int2str(kk) '.mat'];
    load(filename)
        
    fM(i,j)=f;
    w1M(i,j)=w(1);
    w2M(i,j)=w(2);
    YM(i,j)=sum(Y==1);
    
    clearvars -except fM w1M w2M YM kk
    kk
end

%% optimal w values
load('coexp_10_5.mat')
wop=zeros(nParams,1);
lx=find(fM==min(fM(:)));
[rx,cx]=ind2sub(size(fM),lx);
wop(1) = w1M(rx,cx);
wop(2) = w2M(rx,cx);

% update potentials
[nodePot,edgePot] = UGM_CRF_makePotentials(wop,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
% decoding a CRF model
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
Y = int32(Y');

sum(Y==1)


noder = textread('data_hq/trueg.txt','%s');

subT=zeros(1000,1);
for n = 1:1000
        [~,j] = ismember(noder(n),genes);
        subT(n,1)=j;
end

a=find(Y==1);
length(intersect(a,subT(:,1)))
