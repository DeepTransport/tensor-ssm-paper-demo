function [mu, L] = stats_weighted(sams, w)
% computed the mean and approximate range of weighted samples
% sams are multidimensional
% the output is to L(compute sams-mu)
w = w/sum(w); % assure w is normalized
mu = sum(sams.* repmat(w, size(sams, 1),1), 2);

L = zeros(size(sams, 1));
for k = 1:size(sams, 1)
    L(k, k) = sqrt(sum((sams(k, :) - mu(k)).^2.*w));
end
end