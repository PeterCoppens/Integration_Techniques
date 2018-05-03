function [ U ] = nulrules( x )
    n = length(x); % n+1 in fact
    U = zeros(n);
    vandermonde = flipud(vander(x)');
    for m = 1:n-1
        vandermonde = vandermonde(1:end-1,:); 
        nullSpace = null(vandermonde);
        nullVector = nullSpace(:,1);
        for i = 1:m-1
            factor = sum(nullVector.*U(:,i));
            nullVector = nullVector - factor*U(:,i);
        end
        U(:,m) = nullVector;
        U(:,m) = U(:,m)/norm(U(:,m));
    end
end

