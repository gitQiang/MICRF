noder = textread('data_hq/trueg.txt','%s');

subT=zeros(1000,1);
for n = 1:1000
        [~,j] = ismember(noder(n),genes);
       
        subT(n,1)=j;
end

a=find(Y==1);
length(intersect(a,subT(:,1)))

netfile='data_hq/brainspan_net_top5.txt';
beness='data_hq/Betweenness_edge_coexp5.txt';
nodef='data_hq/hotnet_inputmetap.txt';
[w f Y]=nenetLearn(netfile,beness,nodef);

netfile='data_hq/StringNew_HPRD_mnet.txt';
beness='data_hq/Betweenness_edge_HPRD.txt';
nodef='data_hq/hotnet_inputmetap.txt';
[w f Y]=nenetLearn(netfile,beness,nodef);