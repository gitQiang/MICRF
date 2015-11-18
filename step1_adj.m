function adj = step1_adj(netfile)

% read the network edges and node
[node1,node2,weights] = textread(netfile,'%s%s%f','delimiter','\t');
genes = union(node1,node2);
nNodes = length(genes);
% Make Adjacency Matrix 
adj = zeros(nNodes); % Symmetric {0,1} matrix containing edges
for r = 1:length(node1)
    [~,i] = ismember(node1(r),genes);
    [~,j] = ismember(node2(r),genes);
    adj(i,j)=1;
    adj(j,i)=1;
    adj(i,i)=0;
end

end
