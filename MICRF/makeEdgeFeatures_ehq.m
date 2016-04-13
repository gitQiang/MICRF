function [Xedge] = makeEdgeFeatures_ehq(X,edgeEnds,edgePot)

[nInstances,nFeatures,nNodes] = size(X);
nEdges = size(edgeEnds,1);

nEdgeFeatures=4;

% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstances,nEdgeFeatures,nEdges);
for i = 1:nInstances
    for e = 1:nEdges
        Xedge(i,1,e) = log(edgePot(1,1,e));
        Xedge(i,2,e) = log(edgePot(1,2,e));
        Xedge(i,3,e) = log(edgePot(2,1,e));
        Xedge(i,4,e) = log(edgePot(2,2,e));
    end
end