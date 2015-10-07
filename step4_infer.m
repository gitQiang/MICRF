function [Y,nps,eps,logZ]=step4_infer(feaf,wop)

load(feaf)
addpath(genpath(pwd))

nodePot0 = nodePot;
edgePot0 = edgePot;

% node and edge features
nNodes = length(genes);
nInstance = 1;
% Make node features
nNodeFeatures = 2;
Xnode = zeros(nInstance, nNodeFeatures, nNodes);
Xnode(1,1,:) = nodePot0(:,1);
Xnode(1,2,:) = nodePot0(:,2);

% Make edge features
nEdges = size(edgeStruct.edgeEnds,1);
nEdgeFeatures=4;
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

%% decoding and inferring    
% make potentials
w=wop;
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);

% decoding a CRF model
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
Y = int32(Y');
      
% infer conditional probability
[nps,eps,logZ] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);

end
