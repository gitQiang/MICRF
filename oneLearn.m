function oneLearn(netfile,beness,nodef,adj,kk,outputfile)

fprintf('%d\n',kk);

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
%F=tiedrank(be)/length(be);
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

% Make edge features
n1 = edgeStruct.edgeEnds(:,1);
n2 = edgeStruct.edgeEnds(:,2);
edges = strcat(cellstr(genes(n1)),'_',cellstr(genes(n2)));
[~,ind1] = ismember(edges,edgesb1);
[~,ind2] = ismember(edges,edgesb2);
sube = ind1 + ind2;
% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstance,nEdgeFeatures,nEdges);
Xedge(1,1,:) = reshape(1+sqrt(Xnode(1,1,n1) .* Xnode(1,1,n2)),[],1).*F(sube);
Xedge(1,2,:) = F(sube);
Xedge(1,3,:) = F(sube);
Xedge(1,4,:) = reshape(1+sqrt(Xnode(1,2,n1) .* Xnode(1,2,n2)),[],1).*F(sube);
[nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct);
nParams = max([nodeMap(:);edgeMap(:)]);

%% training 
% run with different initial w values
w = ones(nParams,1);
f = 100000;

i = floor((kk-1)/10) + 1;
j = mod(kk,10);
if j==0 
    j = 10; 
end

w(1)=i/10;
w(2)=j/10;
% initial parameters
w0 = w;
f0 = f;
flag = 0;
iter = 0;
inferFunc = @UGM_Infer_LBP; %inferFunc = @UGM_Infer_TRBP; %inferFunc = @UGM_Infer_MeanField; %inferFunc = @UGM_Infer_Junction;  %inferFunc = @UGM_Infer_Block_MF;   %inferFunc = @UGM_Infer_Conditional; 

maxiter=100;
maxFunEvals = 20;
options = [];
options.maxFunEvals = maxFunEvals;
options.progTol = 1e-3;
options.optTol = 1e-3;

while flag==0
    % update potentials
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
    % decoding a CRF model
    Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    Y = int32(Y');

    % training a CRF model
    lambda = ones(size(w)); %lambda(2) = 10;
    regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
    [w,f]= minFunc(regFunObj,w,options);
    
    % fix bugs for illegal direction
    if isnan(f)==1 
        fprintf('%d Here\n',iter);
        w(1)=w(1)+6;
        w(2)=w(2)+26;
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
        Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
        Y = int32(Y');
        lambda = ones(size(w)); %lambda(2) = 10;
        regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
        [w,f]= minFunc(regFunObj,w,options);
    end
    
    if isnan(f)
       w = w0;
       f = f0;
       flag = 2;
       break;
    end
    
    if norm(w-w0,1) <= 1e-5 && abs(f-f0) <= 1e-5 
        flag = 1;
    else
        f0 = f;
        w0 = w;
    end

    iter = iter + 1;
    fprintf('%d\n',iter);
end

fileID = fopen(outputfile,'w');
fprintf(fileID,'%d\t%d\t%d\t%f\t%f\t%f\t%d\n',kk,i,j,w(1),w(2),f,flag);
fclose(fileID);

end
