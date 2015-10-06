function step4_infer()

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


end
