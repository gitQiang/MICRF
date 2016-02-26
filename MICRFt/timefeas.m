function [nodet,edget,t] = timefeas(texp,expG,genes,edgeStruct)
    t = length(texp);
    n1 = edgeStruct.edgeEnds(:,1);
    n2 = edgeStruct.edgeEnds(:,2);
    
    nEdges = edgeStruct.nEdges;
    nodet= 1e-4 * ones(length(genes),t);
    edget= 1e-4 * ones(nEdges,t);
    
    inG=intersect(expG,genes);
    for i=1:t
        rg=mean(texp{i},2);
        rg=tiedrank(rg)/length(expG);
        %rg=(rg-min(rg)+0.01)/(max(rg)-min(rg));
        [~,ind1] = ismember(inG,genes);
        [~,ind2] = ismember(inG,expG);
        nodet(ind1,i)=rg(ind2);
        
        
        subs1=ismember(n1,ind1);
        subs2=ismember(n2,ind1);
        subs = subs1 & subs2;
        
        [~,xsubs]=ismember(genes(n1(subs)),expG);
        [~,ysubs]=ismember(genes(n2(subs)),expG);
        
        X=texp{i}(xsubs,:);
        Y=texp{i}(ysubs,:);
        
        cors=zeros(size(X,1),1);
        for j=1:size(X,1)
            cors(j) = abs(corr2(X(j,:),Y(j,:)));
        end
        
        edget(subs,i)=cors;
    end
        
end
