% first check for the infer results based on top 1000 TADA genes

noder = textread('data_hq/trueg.txt','%s');

subT=zeros(1000,1);
for n = 1:1000
        [~,j] = ismember(noder(n),genes);
        subT(n,1)=j;
end

a=find(Y==1);
length(intersect(a,subT(1:100,1)))
length(intersect(a,subT(1:200,1)))
length(intersect(a,subT(1:500,1)))


%% steps:

% step 1: generate the adj matrix for each network
adj = step1_adj(netfile);
%save  a.mat adj


% step 2: make the node feature and edge feature
step2_feature(netfile,beness,nodef, filename)



% step 3: optimal the parameters

% qsub.m
% getonePara.m



% step 4: give the infer ande decoding results based on optimal parameters



% 27
netfile='data_hq/brainspan_net_top5.txt';
beness='data_hq/Betweenness_edge_coexp5.txt';
nodef='data_hq/hotnet_inputmetap.txt';
[w f Y]=nenetLearn(netfile,beness,nodef);

% 20
netfile='data_hq/StringNew_HPRD_mnet.txt';
beness='data_hq/Betweenness_edge_HPRD.txt';
nodef='data_hq/hotnet_inputmetap.txt';
[w f Y]=nenetLearn(netfile,beness,nodef);


% 3
netfile='data_hq/STRINGnetmap.txt';
beness='data_hq/Betweenness_edge_STRING.txt';
nodef='data_hq/CRF_inputmetap.txt';
[w f Y]=nenetLearn(netfile,beness,nodef);


% 6
netfile='data_hq/iRefIndexm.txt';
beness='data_hq/Betweenness_edge_iRef.txt';
nodef='data_hq/hotnet_inputmetap.txt';

% 7
netfile='data_hq/ASDexp_net_top5.txt';
beness='data_hq/Betweenness_edge_coexp.txt';
nodef='data_hq/hotnet_inputmetap.txt';

% 21
netfile='data_hq/ASDexp_net.txt';
beness='data_hq/Betweenness_edge_coexp1.txt';
nodef='data_hq/hotnet_inputmetap.txt';

% 26
netfile='data_hq/brainspan_net_cor.txt';
beness='data_hq/Betweenness_edge_corr1.txt';
nodef='data_hq/hotnet_inputmetap.txt';

% 28
netfile='data_hq/brainspan_net_top7.txt';
beness='data_hq/Betweenness_edge_coexp7.txt';
nodef='data_hq/hotnet_inputmetap.txt';