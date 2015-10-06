    %% important edge potential definition
    edgePot(1,1,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) = abs(nodePot(n1,1) - nodePot(n2,1))*F(sube);
    edgePot(2,1,e) = abs(nodePot(n1,2) - nodePot(n2,2))*F(sube);
    edgePot(2,2,e) = (1-abs(nodePot(n1,2)-nodePot(n2,2)))*F(sube);
    
        %% important edge potential definition
    edgePot(1,1,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) = (1 - abs(nodePot(n1,2) - nodePot(n2,1)))*F(sube);
    edgePot(2,2,e) = (1 - abs(nodePot(n1,2) - nodePot(n2,2)))*F(sube);

    
    edgePot(1,1,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) = abs(nodePot(n1,1) - nodePot(n2,1))*F(sube);
    edgePot(2,1,e) = abs(nodePot(n1,2) - nodePot(n2,2))*F(sube);
    edgePot(2,2,e) = (1-abs(nodePot(n1,2)-nodePot(n2,2)))*F(sube);
    
    
    edgePot(1,1,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) = (1 - abs(nodePot(n1,2) - nodePot(n2,1)))*F(sube);
    edgePot(2,2,e) = (1 - abs(nodePot(n1,2) - nodePot(n2,2)))*F(sube);
    
    
    edgePot(1,1,e) = (abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(1,2,e) = (1 - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) = (1 - abs(nodePot(n1,2) - nodePot(n2,1)))*F(sube);
    edgePot(2,2,e) = (abs(nodePot(n1,2) - nodePot(n2,1)))*F(sube);
    
    
    n1 = edgeStruct.edgeEnds(e,1);
	n2 = edgeStruct.edgeEnds(e,2);
    sube = 47529;
    
    edgePot(1,1,e) = (nodePot(n1,1) + nodePot(n2,1)) * (1 - abs(nodePot(n1,1) - nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) = (nodePot(n1,1) + nodePot(n2,2)) * (1 - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) = (nodePot(n1,2) + nodePot(n2,1)) * (1 - abs(nodePot(n1,2) - nodePot(n2,1)))*F(sube);
    edgePot(2,2,e) = (nodePot(n1,2) + nodePot(n2,2)) * (1 - abs(nodePot(n1,2) - nodePot(n2,2)))*F(sube);
    
    edgePot(1,1,e) =  - abs(nodePot(n1,1) - nodePot(n2,1))*F(sube);
    edgePot(1,2,e) =  - abs(nodePot(n1,1) - nodePot(n2,2))*F(sube);
    edgePot(2,1,e) =  - abs(nodePot(n1,2) - nodePot(n2,1))*F(sube);
    edgePot(2,2,e) =  - abs(nodePot(n1,2) - nodePot(n2,2))*F(sube);
    
        edgePot(1,1,e) =   abs(nodePot(n1,1) - nodePot(n2,1))*F(sube);
    edgePot(1,2,e) =   abs(nodePot(n1,1) - nodePot(n2,2))*F(sube);
    edgePot(2,1,e) =   abs(nodePot(n1,2) - nodePot(n2,1))*F(sube);
    edgePot(2,2,e) =   abs(nodePot(n1,2) - nodePot(n2,2))*F(sube);
    
    edgePot(1,1,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) =   0.5 * (abs(nodePot(n1,1) - nodePot(n2,2)) - min(abs(nodePot(n1,1)),abs(nodePot(n2,2))))*F(sube);
    edgePot(2,1,e) =   0.5 * (abs(nodePot(n1,2) - nodePot(n2,1)) - min(abs(nodePot(n1,2)),abs(nodePot(n2,1))))*F(sube);
    edgePot(2,2,e) =   0.5 * (abs(nodePot(n1,2)) + abs(nodePot(n2,2)))*F(sube);
    
    edgePot(1,1,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)))*F(sube);
edgePot(1,2,e) =   0.5 * (abs(nodePot(n1,1) - nodePot(n2,1)) - min(abs(nodePot(n1,1)),abs(nodePot(n2,1))))*F(sube);
edgePot(2,1,e) =   edgePot(1,2,e);
edgePot(2,2,e) =   edgePot(1,1,e);

edgePot(1,1,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)) - abs(nodePot(n1,1) - nodePot(n2,2)))*F(sube);
    edgePot(2,1,e) =   edgePot(1,2,e);
    edgePot(2,2,e) =   edgePot(1,1,e);
    
    
        edgePot(1,1,e) =   0.5 * (abs(nodePot(n1,1)) + abs(nodePot(n2,1)))*F(sube);
    edgePot(1,2,e) =   0.5 * (abs(nodePot(n1,1) - nodePot(n2,2)) - min(abs(nodePot(n1,1)),abs(nodePot(n2,2))))*F(sube);
    edgePot(2,1,e) =   0.5 * (abs(nodePot(n1,2) - nodePot(n2,1)) - min(abs(nodePot(n1,2)),abs(nodePot(n2,1))))*F(sube);
    edgePot(2,2,e) =   0.5 * (abs(nodePot(n1,2)) + abs(nodePot(n2,2)))*F(sube);

