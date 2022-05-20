function W_out = nn_graph(n,r) %Outputs graph structure as matrix W_out, and theoretical bifurcation values bif_min, bif_max
    %Create graph structure in matrix form of size n x n
     
    %r = range of interaction scaled between 0 and 1
    W_out = zeros(n,n);
    
    for i = 2:n
        for j = 1:i-1
            if min(abs(i-j)/n,1-abs(i-j)/n)<r
                W_out(i,j) = 1;
            else
                W_out(i,j) = 0;
            end
        end
    end
    W_out = W_out+W_out';
end