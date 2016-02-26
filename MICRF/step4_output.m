function step4_output(genes,Y,nps,Xnode,outputfile,idm)
% step 4: output files
if idm == 1
	fcon =  fopen('/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/Fmap0121.txt','r');
	C = textscan(fcon,'%s%s','delimiter','\t');
	id1=C{1};
	id2=C{2};
	fclose(fcon);
    [Lia,j]=ismember(genes,id1);
    subs = find(Lia > 0);
    sub2 = find(~strcmp(id2(j(subs)),'NA'));
    genes(subs(sub2))=id2(j(subs(sub2)));
end

nodescore = [nps(:,1),reshape(Xnode(1,1,:),[],1)];
% nps with size nNodes * nstate
[~,Ind]=sortrows(-1*nodescore); % 1 as risk, 2 as non-risk
fileID = fopen(outputfile,'w');
for i = 1:length(genes)
fprintf(fileID,'%s\t%12.8f\t%12.8f\t%12.8f\t%d\n',genes{Ind(i)},nps(Ind(i),1),nps(Ind(i),2),Xnode(1,1,Ind(i)),Y(Ind(i)));
end
fclose(fileID);
end