function [w,f,Y]=oneLearn(netfile,beness,nodef,adjf)

load(adjf);
%% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
nStates = ones(1,nNodes) * 2;
maxState = max(nStates);

% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates); 
nEdges = edgeStruct.nEdges;

% read TADA score
fcon =  fopen(nodef,'r');
C = textscan(fcon,'%s%f','delimiter','\t');
nodes= C{1};
score1=C{2};
fclose(fcon);
node_s = intersect(nodes,genes);
% read edge betweenness
[enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
F=tiedrank(be)/length(be);

edgesb1 = strcat(enode1,'_',enode2);
edgesb2 = strcat(enode2,'_',enode1);

% make pseudo potential functions
% 1 risk state, 2 non-risk state
Pri=zeros(2,1);
Pri(1)=0.06;
Pri(2)=0.94;
% Make (non-negative) potential of each node taking each state
nodePot = ones(nNodes,maxState); %!!!!!!!!! zeros and ones for node information
nodePot(:,1) = Pri(1)*Pri(2);
nodePot(:,2) = Pri(2)*Pri(1);
[~,j] = ismember(node_s,genes);
[~,i] = ismember(node_s,nodes);        
nodePot(j,1) = score1(i)*Pri(2);
nodePot(j,2) = (1 - score1(i))*Pri(1);

% Make (non-negative) potential of each edge taking each state combination
edgePot = ones(maxState,maxState,edgeStruct.nEdges);
n1 = edgeStruct.edgeEnds(:,1);
n2 = edgeStruct.edgeEnds(:,2);
edges = strcat(cellstr(genes(n1)),'_',cellstr(genes(n2)));
[~,ind1] = ismember(edges,edgesb1);
[~,ind2] = ismember(edges,edgesb2);
sube = ind1 + ind2;
%% important edge potential definition
edgePot(1,1,:) =   (1+sqrt(nodePot(n1,1) .* nodePot(n2,1))).*F(sube);
edgePot(1,2,:) =   F(sube);
edgePot(2,1,:) =   F(sube);
edgePot(2,2,:) =   (1+sqrt(nodePot(n1,2) .* nodePot(n2,2))).*F(sube);

clear adj
save BrainCor_sim.mat 
