function yval = TT_eval(TT_new, poly, Sigma_update, x)
%   eval the squared convolved TT
%   TT_new is the expanded cores from TT_expand
%   poly is the new basis functions for TT_new
%   x is the input values of size d*n
%   yval is of size 1*n

% m = size(Sigma, 1);
Sigma = Sigma_update(x);
dm = size(x, 1);
nx = size(x, 2);
ind = (any(x <  poly.oned_domains{1}.bound(1))) | ( any(x > poly.oned_domains{1}.bound(2)));

order = poly.oneds{1}.order;
bound = poly.oned_domains{1}.bound(2);
freq2 = [0, (1:order)/bound, (1:order)/bound, (order+1)/bound].^2;
% freq2 = [0, (1:order)/bound, (1:order)/bound].^2;

% modifier(dm-m+1:dm, :,:) = exp(-0.5*pi^2*repmat(diag(Sigma), 1, nx).* reshape(freq2, 1,1,[]));
modifier = exp(-0.5*pi^2* Sigma.* reshape(freq2, 1,1,[]));

phi_x = eval_basis(poly.oneds{1}, x(:)/bound)/bound/2.*reshape(modifier, dm*nx, []);
% phi_x = eval_basis(poly.oneds{1}, x(:)/bound)/bound/2;

x_all = reshape(phi_x * poly.oneds{1}.node2basis, dm, nx, []);
% l = size(x_all, 3);

yval = zeros(1, nx);
for k = 1:nx
    if ind(k) == 1
        yval(k) = 0;
        continue
    end
    TT_new2 = TT_new;
    for k1 = 1:dm
        TT_new2{k1} = TT_new2{k1}.*reshape(x_all(k1, k, :), 1, []);
    end
    yval(k) = ttdot(TT_new, TT_new2);
end

yval = max(yval, 0); % eliminate the roundoff error
% old version
% x_all = reshape(sqrt(phi_x * poly.oneds{1}.node2basis), dm, nx, []);
% l = size(x_all, 3);

% the second dimension is nx, when computing cores, pick x_all(:,k,:)'

% yval = cellfun(@(x) compute(TT_new, x), mat2cell(permute(x_all, [3,1,2]), l, dm, ones(1, nx)));
% yval = reshape(yval, 1, numel(yval));
% yval(1, ind) = 0;
end

