function q = computeq(sams, p)
% sams of size N * T
% compute the p-quantiles for each column of sams, p is a vector
% q is a vector of same size as p

[N, T] = size(sams);
nq = length(p);
q = zeros(nq, T);

for t = 1:T
    temp = sort(sams(:, t), 'ascend');
    for k = 1:nq
        q(k, t) = temp(round(N*p(k)));
    end
end
end