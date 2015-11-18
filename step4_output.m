 function step4_output(netfile,beness,nodef,netflag,onefile,outputfile,w0,idm)

 %% step 1: network adjacency matrix 
 %step1_adj(netfile)
 adjfiles=cell(6,1);
 adjfiles{1}='adj_STRING.mat'; 
 adjfiles{2}='adj_iRef.mat'; 
 adjfiles{3}='adj_HPRD.mat'; 
 adjfiles{4}='adj_braincor.mat'; 
 adjfiles{5}='adj.mat'; 
 adjfiles{6}='adj_Co_PrePPI.mat';  
 adjf=adjfiles{netflag};

 %% step 2: compute initial edge and node features
 genes=step2_feature(netfile,beness,nodef,adjf,onefile)
 
 %% step 3: MICRF training, decoding and inferring
 [Y,nps]=step3_MICRF(onefile,w0)
 
 %% step 4: output files
 if idm == 1
	fcon =  fopen('/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/Fmap0121.txt','r');
	C = textscan(fcon,'%s%s','delimiter','\t');
	id1=C{1};
	id2=C{2};
	fclose(fcon);
       
       for i = 1:length(genes)
       	[Lia,j]=ismember(genes(i),id1)
       	if Lia
       		if ~strcmp(id2(j),'NA')  genes(i)=id2(j); end
       	end
       end         
 end
 [~,Ind]=sort(nps);
fileID = fopen(outputfile,'w');
for i = 1:length(genes)
fprintf(fileID,'%s\t%f\t%d\n',genes{Ind(i)},nps(Ind(i)),Y(Ind(i)));
end
fclose(fileID);

 end