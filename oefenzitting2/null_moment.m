function U = null_moment(q, K)
    % calculate orthogonal equally strong nullrules for
    % the quadrature rule defined by the vector of weights q.
    k = length(q)-1;
    x = -1 + (2/k)*(0:k);
    
    V = flipud(vander(x)');
    U = zeros(k+1, k);
    for m = 1:K
        % Calculate next nullrule
        V = V(1:end-1,:); 
        NS = null(V);
        nv = zeros(size(NS(:, 1)));
        for l = 1:size(NS, 2)
            if norm(U(:, 1:m-1)*(U(:, 1:m-1)\NS(:, l)) - NS(:, l)) > 10^-3
                nv = NS(:, l);
                break;
            end
        end
        
        % Orthogonalize
        for i = 1:m-1
            r = sum(nv.*U(:,i));
            nv = nv - r*U(:,i);
        end
        
        % Make equally strong
        U(:,m) = norm(q)*nv/norm(nv);
    end
    U = U'; % null vectors are stored on the rows
end

