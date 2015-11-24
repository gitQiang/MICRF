function oneLearn(kk,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,outputfile)

fprintf('%d\n',kk);

%% training 
nParams = max([nodeMap(:);edgeMap(:)]);
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

while (flag==0 && iter <= maxiter)
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
