function demoMICRF(netfile,beness,nodef,adjf,feaf)

%% steps:

% step 1: generate the adj matrix for each network
adj = step1_adj(netfile);
% save(adjf,adj)

% step 2: make the node feature and edge feature
step2_feature(netfile,beness,nodef,adjf,feaf)

% step 3: optimal the parameters

% qsub.m
% getonePara.m

% step 4: give the infer ande decoding results based on optimal parameters
[Y,nps,eps,logZ]=step4_infer(feaf,wop)

end