function adj = step1_adjt(node1,node2,genes)
    nNodes = length(genes);
    % Make Adjacency Matrix 
    adj = zeros(nNodes); % Symmetric {0,1} matrix containing edges
    [~,i] = ismember(node1,genes);
    [~,j] = ismember(node2,genes);
    sub1 = nNodes*(j-1) + i;
    sub2 = nNodes*(i-1) + j; 
    adj(sub1)=1;
    adj(sub2)=1;  
    adj(1:(nNodes+1):end)=0;
end
