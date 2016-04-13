function [Xedge] = makeEdgeFeatures_hq(X,edgeEnds,edgePot)

[nInstances,nFeatures,nNodes] = size(X);
nEdges = size(edgeEnds,1);

nEdgeFeatures=1;

% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstances,nEdgeFeatures,nEdges);
for i = 1:nInstances
    for e = 1:nEdges
        %n1 = edgeEnds(e,1);
        %n2 = edgeEnds(e,2);
        %Xedge(i,nEdgeFeatures,e) = log(edgePot(Y(n1),Y(n2)));
        Xedge(i,nEdgeFeatures,e) = log(edgePot(1,1,e));
    end
end