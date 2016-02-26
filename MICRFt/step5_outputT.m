function [coreg,speg,keyg]=step5_outputT(outT,texp,expG,t,outputfile,flagp)
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
    
    % generate key genes from time 2:t: based on rule 1 or rule 2
    % rule 1: mean expressed value decreasing or increasing at least half
    % rule 2: remove redundant genes: 98% genes have the same cor > 0.7
    ro=0.2;
    rhoc=0.7;
    reduc=0.9;
    keyg=cell(t-1,1);
    for i=2:t
        tmpg0=union(tgs{i},tgs{i-1});
        
        % gene expressed values changing
        [~,subs]=ismember(tmpg0,expG);
        tmpg=tmpg0(subs>0);
        subs(subs==0)=[];
        tmp1=mean(texp{i-1}(subs,:),2);
        tmp2=mean(texp{i}(subs,:),2);
        subs=(tmp1./tmp2 <= 1-ro) | (tmp1./tmp2 >= 1+ro);
        if all(subs==0)
            ekg=[];
        else
            ekg=tmpg{subs>0};
            ekg=cellstr(ekg);
        end
        
        % gene rank values changing
        [~,rnd1]=ismember(tmpg0,outT{i-1}(:,1));
        [~,rnd2]=ismember(tmpg0,outT{i}(:,1));
        subs=(rnd1./rnd2 <= 1-ro);
        if all(subs==0)
            rkg=[];
        else
            rkg=tmpg0{subs>0};
            rkg=cellstr(rkg);
        end
        
        %onep=union(ekg,rkg);
        onep=setdiff(tgs{i},tgs{i-1});
        if isempty(onep)
            keyg{i-1}=setdiff(tgs{i},tgs{i-1});
        elseif length(onep)==1
            keyg{i-1}=onep;
        else
            % delete redundant co-exp genes
            % order onep
            np=length(onep);
            [~,ind]=ismember(onep,outT{i}(:,1));
            onep=outT{i}(sort(ind),1);
            
            X=texp{i}';
            coEg=cell(np,1);
            for ii=1:np
                [~,ind]=ismember(onep(ii),expG);
                rho=corr(X,texp{i}(ind,:)');
                coEg{ii}=expG(rho>rhoc);
            end
            
            onepNew=onep(1);
            for ii=2:np
                flag=0;
                for jj=1:(ii-1)
                    tmp=intersect(coEg{ii},coEg{jj});
                    if length(tmp)/length(coEg{ii}) >= reduc
                        flag=1;
                    end
                end
                if flag==0
                    onepNew=union(onepNew,onep(ii));
                end
            end
            keyg{i-1}=onepNew;
        end
    end
    
    if strcmp(outputfile,'') == 0
        outputfile1=[outputfile,'Summary.txt'];
        fileID = fopen(outputfile1,'w');
        fprintf(fileID,'# core risk genes across all time points.\n');
        for i = 1:length(coreg)
            fprintf(fileID,'%s\n',coreg{i});
        end
        
        fprintf(fileID,'# Time point specific risk genes. \n');
        for i = 1:t
            fprintf(fileID,'# Specific risk genes at time %d .\n',i);
            for j=1:length(speg{i})
                fprintf(fileID,'%s\n',speg{i}{j});
            end
        end
        
        fprintf(fileID,'# Time points transformed risk genes. \n');
        for i = 1:(t-1)
            fprintf(fileID,'# Risk genes transformed from time %d to time %d.\n',i,i+1);
            for j=1:length(keyg{i})
                fprintf(fileID,'%s\n',keyg{i}{j});
            end
        end
        
        fclose(fileID);
    end
    
    if flagp
        %intX=2;
        intT=20;
        X=[];
        Y=[];
        col=[];
        xla=[];
        str=[];
        fsize=[];
        fsize1=[];
        for i=2:t
            [~,xone]=ismember(keyg{i-1},outT{i}(:,1));
            % jitter xone
            %psex=1:intX:(length(xone)*intX);
            %xone=psex(tiedrank(xone));
            
            tmp=mean(texp{i},2);
            %tmpR=tiedrank(tmp);
            [~,ind]=ismember(keyg{i-1},expG);
            %yone=tmpR(ind)/length(expG);
            yone=tmp(ind);            
            
            if isempty(X) xmax=0; else xmax=max(X)+intT; end
            X=[X; reshape(xone+xmax,[],1)-3]; % position with circles
            Y=[Y; reshape(yone,[],1)];
            col=[col; i*ones(length(yone),1)];
            xla=[xla median(reshape(xone+xmax,[],1))];
            str=[str cellstr(keyg{i-1})];
            fsize=[fsize;reshape((max(xone)+1-xone)/max(xone),[],1)];
            fsize1=[fsize1;reshape(xone,[],1)];
        end
        
        hFig = figure();
        set(hFig, 'Visible', 'off')
        scatter(X,Y,15,col);
        for ii=1:size(X,1)
            text(X(ii),Y(ii),str(ii),'FontSize',floor(fsize(ii)*5+4));
        end
        ylim([min(Y)-0.7 max(Y)+0.5])
        xlabel('Time-series Changing')
        ylabel('Gene expression level')
        ax=gca;
        set(ax,'XTick',xla,'XTickLabel',{'1-2','2-3','3-4','4-5','5-6'})
        print(gcf, '-dpdf', 'MICRFtPlotKeys.pdf')
        
        % write plot information
        fileID = fopen('MICRFtPlotKeys.txt','w');
        for ii=1:size(X,1)
            fprintf(fileID,'%12.8f\t%12.8f\t%s\t%f\t%d\n',X(ii),Y(ii),str{ii},fsize1(ii),col(ii));
        end
        fclose(fileID);
        
        
        %plot core genes across time points
        X=[];
        Y=[];
        [~,ind]=ismember(coreg,expG);
        for i=1:t
            [~,xone]=ismember(coreg,outT{i}(:,1));
            tmp=mean(texp{i},2);
            yone=tmp(ind);
            X=[X,xone];
            Y=[Y,yone];  
        end
        
        x=mean(X,2);
        y=mean(Y,2);
        fsize=(max(x)-x+1)/max(x);
        hFig = figure();
        set(hFig, 'Visible', 'off')
        scatter(x,y,15);
        for ii=1:length(x)
            text(x(ii),y(ii),coreg(ii),'FontSize',4);
        end
        ylim([min(y)-0.7 max(y)+0.5])
        xlabel('Mean Ranks of core genes')
        ylabel('Gene expression level')
        ax=gca;
        print(gcf, '-dpdf', 'MICRFtPlotCores.pdf')
        
        % write plot information
        fileID = fopen('MICRFtPlotCores.txt','w');
        for ii=1:length(x)
            fprintf(fileID,'%12.8f\t%12.8f\t%s\n',x(ii),y(ii),coreg{ii});
        end
        fclose(fileID);
        
    end

end