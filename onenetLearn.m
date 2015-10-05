function [w,f,Y]=onenetLearn(netfile,beness,nodef)

%% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
nStates = ones(1,nNodes) * 2;
maxState = max(nStates);

%% Make Adjacency Matrix and EdgeStruct
%nStates: % Number of states that each node can take
% maxState % Maximum number of states that any node can take
% nNodes  % Total number of nodes

% adj = zeros(nNodes); % Symmetric {0,1} matrix containing edges
% for r = 1:length(node1)
%     [~,i] = ismember(node1(r),genes);
%     [~,j] = ismember(node2(r),genes);
%     adj(i,j)=1;
%     adj(j,i)=1;
%     adj(i,i)=0;
% end

% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates); 
nEdges = edgeStruct.nEdges;

%% Make the non-negative node and edge potentials
% read TADA score
[nodes, score1, score2] = textread(nodef,'%s%f%f','delimiter','\t');
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

% Make (non-negative) potential of each node taking each state
nodePot = zeros(nNodes,maxState); %!!!!!!!!! zeros and ones for node information
nodePot(:,1) = log10(0.06);
nodePot(:,2) = -log10(0.06);
for n = 1:length(node_s)
        [~,j] = ismember(node_s(n),genes);
        [~,i] = ismember(node_s(n),nodes);
        
        nodePot(j,1) = score1(i);
        nodePot(j,2) = score2(i);
        %nodePot(j,:) = nodePot(j,:)/max(nodePot(j,:));
end

% amax = max(nodePot(:,1));
% amin = min(nodePot(:,1));
% nodePot(:,1)= (nodePot(:,1)-amin)/(amax-amin);
% 
% amax = max(nodePot(:,2));
% amin = min(nodePot(:,2));
% nodePot(:,2)= (nodePot(:,2)-amin)/(amax-amin);


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
    edgePot(1,1,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)) - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) =   edgePot(1,2,e);
    edgePot(2,2,e) =   edgePot(1,1,e);
end

%% loop until convergence 
nodePot0 = nodePot;
edgePot0 = edgePot;


%% initial node and edge features
addpath(genpath(pwd))
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


%% run with different initial w values

w = ones(nParams,1);
fM = -1 * ones(10,10);
w1M = -1 * ones(10,10);
w2M = -1 * ones(10,10);
for i=1:10
    for j=1:10
        w(1)=i/10;
        w(2)=j/10;
        % initial parameters
        w0 = w;
        flag = 0;
        iter = 0;
        inferFunc = @UGM_Infer_LBP; %inferFunc = @UGM_Infer_TRBP; %inferFunc = @UGM_Infer_MeanField; %inferFunc = @UGM_Infer_Junction;  %inferFunc = @UGM_Infer_Block_MF;   %inferFunc = @UGM_Infer_Conditional; 

        while flag==0
        % update potentials
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
        % decoding a CRF model
        Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
        Y = int32(Y');
        % training a CRF model
        options = struct('Display','iter','MaxIter',100,'TolX',1e-5);
        lambda = 10*ones(size(w)); %lambda(2) = 10;
        regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
        [w,f]= minFunc(regFunObj,w);
        %[w,f] = minFunc(@UGM_CRF_NLL,w,options,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);

        if norm(w-w0,1) <= 1e-5
            flag = 1;
        else
            w0 = w;
        end

        iter = iter + 1;
        fprintf('%d\n',iter);
        end
        
        fM(i,j)=f;
        w1M(i,j)=w(1);
        w2M(i,j)=w(2);
    end
end

%% optimal w values
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


%% conditional probability for each gene

