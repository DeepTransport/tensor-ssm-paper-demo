function v = compute_tensor_prod(cores, x)
%compute the evaluation of squared convolved TT specified by cores at x
%   cores are the squared convolved TT of length d+m
%   x of size l*(d+m)
d = size(cores, 1);

% start from the last core d
[r0, ~] = size(cores{d});
core = cores{d}  .*repmat(x(:, d)', r0, 1);

v = reshape(core*core.', [], 1);

for k = d-1:-1:1
    [r0, ~, r1] = size(cores{k});
    core = cores{k}.*repmat(x(:, k)', r0, 1, r1);
    A = zeros(r0^2,r1^2);
    for k1 = 1:r0^2
        k1a = ceil(k1/r0);
        k1b = k1 - r0*k1a + r0;
        for k2 = 1:r1^2
           k2a = ceil(k2/r1);
           k2b = k2 - r1*k2a + r1;
           A(k1, k2) =  core(k1a, :, k2a)*core(k1b,:,k2b).';
        end
    end
    v = A*v;
    
%     core_new = reshape(reshape(core, r0*l, r1)*L, r0, l*r1);
%     if k ~= 1
%         L = chol(real(core_new*core_new'))';
%     else
%         yval = real(core_new*core_new');
%     end
end

end

