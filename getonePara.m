function wop=getonePara(onefile,dirstr)

fM = -1 * ones(10,10);
w1M = -1 * ones(10,10);
w2M = -1 * ones(10,10);
YM = -1 * ones(10,10);

for kk=1:100

    i = floor((kk-1)/10) + 1;
    j = mod(kk,10);
    if(j==0) j = 10; end
    
    filename= [dirstr int2str(kk) '.mat'];
    load(filename)
        
    fM(i,j)=f;
    w1M(i,j)=w(1);
    w2M(i,j)=w(2);
    YM(i,j)=sum(Y==1);
    
    clearvars -except fM w1M w2M YM kk dirstr onefile
    kk
end

%% optimal w values
addpath(genpath(pwd))
load(onefile)

nodePot0 = nodePot;
edgePot0 = edgePot;
%% initial node and edge features
nNodes = length(genes);
nInstance = 1;
% Make node features
nNodeFeatures = 2;
Xnode = zeros(nInstance, nNodeFeatures, nNodes);
Xnode(1,1,:) = nodePot0(:,1);
Xnode(1,2,:) = nodePot0(:,2);

% Make edge features
%Xedge = makeEdgeFeatures_ehq(Xnode,edgeStruct.edgeEnds,edgePot0);
nEdges = size(edgeStruct.edgeEnds,1);
nEdgeFeatures=4;
% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstance,nEdgeFeatures,nEdges);
for i = 1:nInstance
    for e = 1:nEdges
        Xedge(i,1,e) = edgePot0(1,1,e);
        Xedge(i,2,e) = edgePot0(1,2,e);
        Xedge(i,3,e) = edgePot0(2,1,e);
        Xedge(i,4,e) = edgePot0(2,2,e);
    end
end

[nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct);
nParams = max([nodeMap(:);edgeMap(:)]);

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
length(intersect(a,subT(1:100,1)))
length(intersect(a,subT(1:200,1)))
length(intersect(a,subT(1:500,1)))

end


