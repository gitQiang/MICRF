function out=step4_outputt(genes,Y,nps,Xnode,outputfile)

% step 4: output files
nodescore = [nps(:,1),reshape(Xnode(1,1,:),[],1)];
% nps with size nNodes * nstate
[~,Ind]=sortrows(-1*nodescore); % 1 as risk, 2 as non-risk
ng = length(genes);
out=cell(ng,5);
for i = 1:ng
    out(i,:)={genes{Ind(i)},nps(Ind(i),1),nps(Ind(i),2),Xnode(1,1,Ind(i)),Y(Ind(i))};
end

if strcmp(outputfile,'') == 0
    fileID = fopen(outputfile,'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n','Gene','Risk-state','Non-risk-state','Node-feature','Label');
    for i = 1:ng
        fprintf(fileID,'%s\t%12.8f\t%12.8f\t%12.8f\t%d\n',genes{Ind(i)},nps(Ind(i),1),nps(Ind(i),2),Xnode(1,1,Ind(i)),Y(Ind(i)));
    end
    fclose(fileID);
end

end
