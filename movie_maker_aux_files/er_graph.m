function W_out = er_graph(n,p) %Outputs graph structure as matrix W_out, and theoretical bifurcation values bif_min, bif_max
    %Create graph structure in matrix form of size n x n
     
        %Erdos-Renyi graph
            %probability of forming an edge, 0 < p < 1
            % Prob( rand() < p ) = p 
            W_rand = rand(n,n);
            W_rand(W_rand < p) = 1;  %Form an edge
            W_rand(W_rand ~= 1) = 0;  %No edge

            %Need to make the matrix symmetric; store the upper triangular part
            %and symmetrize the resultant matrix       
            W_upper = triu(W_rand,1);

            W_out = diag(zeros(n,1))+W_upper+W_upper';

end