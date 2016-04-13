function [coreg,speg,out]=step5_outputT(outT,t,outputfile)
    tgs=cell(t,1);
    for i=1:t
        tgs{i}=outT{i}(cell2mat(outT{i}(:,2))==1,1);
    end

    coreg=tgs{1};
    for i=2:t
        coreg=intersect(coreg,tgs{i});
    end
    
    speg=cell(t,1);
    for i=1:t
        speg{i}=setdiff(tgs{i},coreg);
    end
    
    % votes from different time-points
    genes = outT{1}(:,1);
    ng = length(genes);
    indM = zeros(ng,t);
    scoreM = zeros(ng,4);
    for i=1:t
        [~,indM(:,i)]=ismember(genes,outT{i}(:,1));
        for j=1:4
            aa=outT{i}(indM(:,i),j+1);
            aa=cell2mat(aa);
            aa=double(aa);
            scoreM(:,j)=scoreM(:,j)+aa;
        end
    end
    sumR = sum(indM,2);
    [~,Ind]=sort(sumR);  
    out=cell(ng,6);
    out(:,1)=genes(Ind);
    out(:,2)=num2cell(sumR(Ind)/t);
    for i=1:4
        out(:,i+2)=num2cell(scoreM(Ind,i)/t);
    end
    
    
    if strcmp(outputfile,'') == 0
    fileID = fopen(outputfile,'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\n','Gene','Time-Sum-Rank','Risk-state','Non-risk-state','Node-feature','Label');
    for i = 1:ng
        fprintf(fileID,'%s\t%f\t%12.8f\t%12.8f\t%12.8f\t%f\n',out{i,1},out{i,2},out{i,3},out{i,4},out{i,5},out{i,6});
    end
    fclose(fileID);
    end

end