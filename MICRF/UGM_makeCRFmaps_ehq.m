function [nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct)
% Assumes that all nodes have the same number of states

nNodes = size(Xnode,3);
nEdges = edgeStruct.nEdges;
nStates = edgeStruct.nStates;
maxState = max(nStates);

UGM_assert(min(nStates)==maxState,'UGM_makeCRFMaps assumes that all nodes must have the same number of states');
nStates = nStates(1);

nNodeFeatures = size(Xnode,2);
nEdgeFeatures = size(Xedge,2);

nodeMap = zeros(nNodes,nStates,nNodeFeatures,'int32');

fNum=1;
nodeMap(:,1,1) = fNum;
nodeMap(:,1,2) = 0;
nodeMap(:,2,1) = 0;
nodeMap(:,2,2) = fNum;

fNum=fNum+1;

nNodeParams = max(nodeMap(:));

edgeMap = zeros(nStates,nStates,nEdges,nEdgeFeatures,'int32');
edgeMap(1,1,:,1) = fNum;
edgeMap(1,2,:,2) = fNum;
edgeMap(2,1,:,3) = fNum;
edgeMap(2,2,:,4) = fNum;

end