function wop0=getOne()

%% different for simulation sets
nSim=100;
outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/all/MICRFs_';
wop0= zeros(2,6);

w1V=zeros(600,1);
w2V=zeros(600,1);
fV =zeros(600,1);

for kk = 1:600
    netj=mod(kk,nSim);
    if netj == 0 
        netj=nSim;
    end
    netflag=floor((kk-1)/nSim)+1;
    outputfile=[outputstr,netflag,'_',netj,'.txt'];

    fcon =  fopen(outputfile,'r');
    C = textscan(fcon,'%d\t%d\t%d\t%f\t%f\t%f\t%d\n','delimiter','\t');
    fclose(fcon);

    w1V(kk) = C{4};
    w2V(kk) = C{5};
    fV(kk)= C(6);
end

for i = 1:6
    f100 = fV(((i-1)*100+1):(i*100));
    sub0 = find(f100==min(f100));
    sub = (i-1)*100 + sub0;
    wop0(1,i)=w1V(sub(1));
    wop0(2,i)=w2V(sub(1));
end

end

