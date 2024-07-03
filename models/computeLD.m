function [mu, L, scaleL] = computeLD(sams, varargin)
% computed the mean and approximate range of weighted samples
% sams are multidimensional
% the output is to L(compute sams-mu)

d = size(sams, 1);
if isempty(varargin)
    w = ones(1, size(sams, 2));
else
    w = varargin{1};
end

q = 0.01;

w = w/sum(w); % assure w is normalized
mu = sum(sams.* repmat(w, d, 1), 2);

cs = sams - mu; % centered samples
Sigma = (cs.*w)*cs';
% L = chol(Sigma + 1e-5* eye(d))';
L = diag(sqrt(diag(Sigma)));


if ESS(w) > 1000
    ss = L\cs; % standarized samples
    % find the stretcher based on the quantiles of ss

    scaleL = zeros(d);
    for k = 1: d
        [sams_order, I] = sort(ss(k, :));
        w_order = w(I);
        b1 = sams_order(cumsum(w_order) > q);
        b2 = sams_order(cumsum(w_order) > 1-q);
        scaleL(k, k) = b2(1) - b1(1);
    end
    scaleL = -scaleL/norminv(q)/2;
else
    scaleL = eye(d);
end
L = L*scaleL;

% disp(scaleL)


% L = zeros(size(sams, 1));
% for k = 1:size(sams, 1)
%     L(k, k) = sqrt(sum((sams(k, :) - mu(k)).^2.*w));
% end
end