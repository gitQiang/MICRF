function [genes,Xnode,Xedge,nodeMap,edgeMap,edgeStruct]=step2_feature(netfile,beness,nodef,adj)

% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
nStates = ones(1,nNodes) * 2;

% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates); 

% read TADA score
fcon =  fopen(nodef,'r');
C = textscan(fcon,'%s%s','delimiter','\t');
nodes =C{1};
score1=C{2};
fclose(fcon);
% delete NA value in score1
subs=find(~strcmp(score1,'NA'));
nodes=nodes(subs);
score1=score1(subs);
score1=str2double(score1);

node_s = intersect(nodes,genes);
% read edge betweenness
[enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
%F=tiedrank(be)/length(be); %uniform distribution
%cuto = quantile(be,0.75)+1.5*(quantile(be,0.75)-quantile(be,0.25));
%be(be>cuto)=cuto;
%F=be/max(be); % log normal distribution
F=log(be)/max(log(be)); % normal distritbution

edgesb1 = strcat(enode1,'_',enode2);
edgesb2 = strcat(enode2,'_',enode1);

%% initial node and edge features: % 1 risk state, 2 non-risk state
% risk genes prior information
Pri=zeros(2,1); 
Pri(1)=0.06;
Pri(2)=0.94;

nInstance = 1;
nNodes = length(genes);
nEdges = size(edgeStruct.edgeEnds,1);
nEdgeFeatures=4;

% Make node features
nNodeFeatures = 2;
Xnode = zeros(nInstance, nNodeFeatures, nNodes);
Xnode(1,1,:) = Pri(1)*Pri(2);
Xnode(1,2,:) = Pri(2)*Pri(1);
[~,j] = ismember(node_s,genes);
[~,i] = ismember(node_s,nodes);        
Xnode(1,1,j) = score1(i)*Pri(2);
Xnode(1,2,j) = (1 - score1(i))*Pri(1);
Xnode = Xnode / mean(Xnode(:));

% Make edge features
n1 = edgeStruct.edgeEnds(:,1);
n2 = edgeStruct.edgeEnds(:,2);
edges = strcat(cellstr(genes(n1)),'_',cellstr(genes(n2)));
[~,ind1] = ismember(edges,edgesb1);
[~,ind2] = ismember(edges,edgesb2);
sube = ind1 + ind2;
% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstance,nEdgeFeatures,nEdges);
Xedge(1,1,:) = reshape(min(Xnode(1,1,n1), Xnode(1,1,n2)),[],1).*F(sube);
Xedge(1,2,:) = reshape(min(Xnode(1,1,n1), Xnode(1,2,n2)),[],1).*F(sube);
Xedge(1,3,:) = reshape(min(Xnode(1,2,n1), Xnode(1,1,n2)),[],1).*F(sube);
Xedge(1,4,:) = reshape(max(Xnode(1,2,n1), Xnode(1,2,n2)),[],1).*F(sube);
[nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct);

end
