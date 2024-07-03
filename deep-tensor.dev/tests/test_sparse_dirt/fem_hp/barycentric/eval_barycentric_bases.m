function f = eval_barycentric_bases(alpha, A, lambdas)
%eval_barycentric_bases
%
% evaluate basis functions for inputs lambda
%
% Tiangang Cui, 29/Oct/2016

% lambdas:  n x nd
% f:        n x poly_order 

n = size(lambdas,1);
f = zeros(n, size(A,1));
for i = 1:n
    f(i,:) = A*prod(repmat(lambdas(i,:), size(alpha,1), 1).^alpha, 2);
end

end