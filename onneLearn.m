function [w,f,Y]=onneLearn(netfile,beness,nodef)

%% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
nStates = ones(1,nNodes) * 2;
maxState = max(nStates);

% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates); 
nEdges = edgeStruct.nEdges;

%% Make the non-negative node and edge potentials
% read TADA score
[nodes, score1] = textread(nodef,'%s%f','delimiter','\t');
node_s = intersect(nodes,genes);
% read edge betweenness
[enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
%be = be /max(be);
n_e = length(be);
F = zeros(1,n_e);
for i=1:n_e
    F(i) = sum(be <= be(i))/n_e;
end


%% make pseudo potential functions
% 1 risk state, 2 non-risk state

Pri=zeros(2,1);
Pri(1)=0.06;
Pri(2)=0.94;
% Make (non-negative) potential of each node taking each state
nodePot = ones(nNodes,maxState); %!!!!!!!!! zeros and ones for node information
nodePot(:,1) = Pri(1)*Pri(2);
nodePot(:,2) = Pri(2)*Pri(1);

for n = 1:length(node_s)
     [~,j] = ismember(node_s(n),genes);
     [~,i] = ismember(node_s(n),nodes);
        
     nodePot(j,1) = score1(i)*Pri(2);
     nodePot(j,2) = (1 - score1(i))*Pri(1);
end


% Make (non-negative) potential of each edge taking each state combination
edgePot = ones(maxState,maxState,edgeStruct.nEdges);
for e = 1:nEdges
	n1 = edgeStruct.edgeEnds(e,1);
	n2 = edgeStruct.edgeEnds(e,2);
    edgePot(:,:,e) = ones(nStates(n1),nStates(n2));
    
    tmp1 = ismember(enode1,genes(n1));
    tmp2 = ismember(enode2,genes(n2));
    tmp = tmp1 + tmp2;
    
    if max(tmp) < 2 
        tmp1 = ismember(enode1,genes(n2));
        tmp2 = ismember(enode2,genes(n1));
        tmp = tmp1 + tmp2;
    end
    sube = find(tmp==2);
    
    %% important edge potential definition
    
    edgePot(1,1,e) =   (1+sqrt(nodePot(n1,1) * nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) =   F(sube);
    edgePot(2,1,e) =   F(sube);
    edgePot(2,2,e) =   (1+sqrt(nodePot(n1,2) * nodePot(n2,2)))*F(sube);
end

clear adj
save Coexp5_sim.mat 
