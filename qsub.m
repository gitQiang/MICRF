function a = qsub(kk)

addpath(genpath(pwd))
load('coexp_max_1_v4.mat')

%% loop until convergence 
nodePot0 = nodePot;
edgePot0 = edgePot;


%% initial node and edge features
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

        while flag==0
        % update potentials
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
        % decoding a CRF model
        Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
        Y = int32(Y');
        % training a CRF model
        options = struct('Display','iter','MaxIter',100,'TolX',1e-5);
        lambda = ones(size(w)); %lambda(2) = 10;
        regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
        [w,f]= minFunc(regFunObj,w);
        %[w,f] = minFunc(@UGM_CRF_NLL,w,options,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);

        if norm(w-w0,1) <= 1e-5
        %if abs(f-f0) <= 1e-5 
            flag = 1;
        else
            f0 = f;
            w0 = w;
        end

        iter = iter + 1;
        fprintf('%d\n',iter);
        end
        
        
        filename =[int2str(kk) '.mat'];
save(filename) 


