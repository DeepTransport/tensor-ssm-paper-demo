function yval = compute(cores, x)
%compute the evaluation of squared convolved TT specified by cores at x
%   cores are the squared convolved TT of length d+m
%   x of size l*(d+m)
d = size(cores, 1);

% start from the last core d
[r0, ~] = size(cores{d});
core = cores{d}.*repmat(x(:, d).', r0, 1);

L = chol(real(core*core.') + 1e-8*eye(r0)).';
% [Q, R] = qr(core.', "econ");
% Q = Q*diag(sign(diag(R)));
% L1 = (diag(sign(diag(R)))*R)';

% if norm(L1-L) > 1e-4 && norm(L1+L) > 1e-4
%     error
% end
% svd(core)

for k = d-1:-1:1
    [r0, l, r1] = size(cores{k});
    core = cores{k}.*repmat(x(:, k).', r0, 1, r1);
    core_new = reshape(reshape(core, r0*l, r1)*L, r0, l*r1);
    if k ~= 1
        L = chol(real(core_new*core_new.') + 1e-8*eye(r0)).';
%         svd(core_new)
        % [Q, R] = qr(core_new.', "econ");
        % Q = Q*diag(sign(diag(R)));
        % L1 = (diag(sign(diag(R)))*R)';
        % if norm(L1-L) > 1e-4 && norm(L1+L) > 1e-4
        %     error
        % end
    else
        yval = real(core_new*core_new.');
    end
end

end

