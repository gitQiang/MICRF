
%beness='data_hq/Betweenness_edge_coexp5.txt';
%beness='data_hq/Betweenness_edge_corr1.txt';
beness='data_hq/Betweenness_edge_Co_PrePPI.txt';
%beness='data_hq/Betweenness_edge_iRef.txt';
%beness='data_hq/Betweenness_edge_STRING.txt';
%beness='data_hq/Betweenness_edge_HPRD.txt';

figure;
[enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
a=fitdist(be,'Lognormal');
histfit(log(be),floor(sqrt(length(be))),'normal')

title('CoPrePPI network betweenness distribution','FontSize',15,'FontWeight','bold')
xlabel('log scale betweenness','FontSize',15,'FontWeight','bold') % x-axis label
ylabel('Frequency','FontSize',15,'FontWeight','bold') % y-axis label

%===============================================================================
netfile='data_hq/brainspan_net_top5.txt';
beness='data_hq/Betweenness_edge_coexp5.txt';
nodef='data_hq/hotnet_inputmeta.txt';
load('adj')

[genes,Xnode,Xedge,nodeMap,edgeMap,edgeStruct]=step2_feature(netfile,beness,nodef,adj);

figure;
histfit(reshape(Xnode(1,1,:),[],1),floor(sqrt(size(Xnode,3))),'Gamma')
title('Node feature distribution (Coexp)','FontSize',15,'FontWeight','bold')
xlabel('Node feature of risk state','FontSize',15,'FontWeight','bold') % x-axis label
ylabel('Frequency','FontSize',15,'FontWeight','bold') % y-axis label

figure;
histfit(reshape(Xnode(1,2,:),[],1),floor(sqrt(size(Xnode,3))),'Gamma')
title('Node feature distribution (Coexp)','FontSize',15,'FontWeight','bold')
xlabel('Node feature of non-risk state','FontSize',15,'FontWeight','bold') % x-axis label
ylabel('Frequency','FontSize',15,'FontWeight','bold') % y-axis label