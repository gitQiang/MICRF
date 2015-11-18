function [Y,nps,eps,logZ]=step3_MICRF(netfile,beness,nodef,adjf,w0)

load(adjf);

% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
nStates = ones(1,nNodes) * 2;

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
n_e = length(be);
F = zeros(n_e,1);
for i=1:n_e
    F(i) = sum(be <= be(i))/n_e;
end
edgesb1 = strcat(enode1,'_',enode2);
edgesb2 = strcat(enode2,'_',enode1);

% risk genes prior information
Pri=zeros(2,1); 
Pri(1)=0.06;
Pri(2)=0.94;

%% initial node and edge features: % 1 risk state, 2 non-risk state
nNodes = length(genes);
nInstance = 1;
% Make node features
nNodeFeatures = 2;
Xnode = zeros(nInstance, nNodeFeatures, nNodes);
Xnode(1,1,:) = Pri(1)*Pri(2);
Xnode(1,2,:) = Pri(2)*Pri(1);
[~,j] = ismember(node_s,genes);
[~,i] = ismember(node_s,nodes);        
Xnode(1,1,j) = score1(i)*Pri(2);
Xnode(1,2,j) = (1 - score1(i))*Pri(1);


% Make edge features
n1 = edgeStruct.edgeEnds(:,1);
n2 = edgeStruct.edgeEnds(:,2);
edges = strcat(cellstr(genes(n1)),'_',cellstr(genes(n2)));
[~,ind1] = ismember(edges,edgesb1);
[~,ind2] = ismember(edges,edgesb2);
sube = ind1 + ind2;
nEdges = size(edgeStruct.edgeEnds,1);
nEdgeFeatures=4;
% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstance,nEdgeFeatures,nEdges);
Xedge(i,1,:) = (1+sqrt(nodePot(n1,1) .* nodePot(n2,1))).*F(sube);
Xedge(i,2,:) = F(sube);
Xedge(i,3,:) = F(sube);
Xedge(i,4,:) = (1+sqrt(nodePot(n1,2) .* nodePot(n2,2))).*F(sube);
[nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct);


%% training 
% initial parameters
f = 100000;
f0 = f;
w=w0;
flag = 0;
iter = 0;
inferFunc = @UGM_Infer_LBP; 
while flag==0
    % update potentials
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
    % decoding a CRF model
    Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    Y = int32(Y');

    % training a CRF model
    lambda = ones(size(w)); 
    regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
    [w,f]= minFunc(regFunObj,w);
   
    if norm(w-w0,1) <= 1e-5 & norm(f-f0,1) <= 1e-5
        flag = 1;
    else
        f0 = f;
        w0 = w;
    end
    
    iter = iter + 1;
    fprintf('%d\n',iter);
end


%% decoding and inferring    
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct); % decoding a CRF model
Y = int32(Y');
[nps,eps,logZ] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct); % infer conditional probability

end
