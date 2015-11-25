function [Y,nps]=step3_MICRF(Xnode,Xedge,nodeMap,edgeMap,edgeStruct,w0)

%% training 
% initial parameters
f = 100000;
f0 = f;
w=w0;
flag = 0;
iter = 0;
maxiter=100;
inferFunc = @UGM_Infer_LBP; 

maxFunEvals = 20;
options = [];
options.maxFunEvals = maxFunEvals;
options.progTol = 1e-3;
options.optTol = 1e-3;
        
while (flag==0 && iter <= maxiter)
    % update potentials
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
    % decoding a CRF model
    Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    Y = int32(Y');

    % training a CRF model
    lambda = ones(size(w)); 
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
        break;
    end
    if norm(w-w0,1) <= 1e-5 && norm(f-f0,1) <= 1e-5
        flag = 1;
    else
        f0 = f;
        w0 = w;
    end
    
    iter = iter + 1;
    fprintf('%d\n',iter);
end

fprintf('%f\t%f\t%f\n',w(1),w(2),f);
%% decoding and inferring    
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct); % decoding a CRF model
Y = int32(Y');
nps = UGM_Infer_LBP(nodePot,edgePot,edgeStruct); % infer conditional probability
if any(isnan(nps))
    % nps = UGM_Infer_TRBP(nodePot,edgePot,edgeStruct); % trying TRBP
    nps(isnan(nps)) = Xnode(isnan(nps)); 
end

end
