function [cores_new, poly_new] = TT_expand(sirt)
%   Computes the cores of TT after expansion
%   sirt is the original tensor train for square root
%   
%   Output cores are the cores containing evaluations(not coefficients!) for the sqaure of TT
%   Length of TT is d+m
%   the degree of freedom is doubled roughly

% m = size(Sigma, 1);
cores_old = sirt.approx.data.cores;
dm = size(cores_old, 1);
polys = sirt.approx.base.oneds;
% bound = sirt.approx.oned_domains{1}.bound(2);
order = polys{1}.order;


poly_new = ApproxBases(Fourier(2*order), sirt.approx.base.oned_domains{1}, dm);

nodes = poly_new.oneds{1}.nodes;
fun_sqr = @(x) eval_radon(polys{1}, x, nodes); % fun input to cellfun
% x is the coeff stored in the old cores, column vector of 2order+1
% output is column vector of length 4order+1

cores_new = cell(dm, 1);

for k = 1:dm % all dimensions, compute square
    core = cores_old{k};
    [r0, n, r1] = size(core);
    cell_core = mat2cell(reshape(permute(core, [2,1,3]), n, []), n, ones(1, r0*r1));
    cell_core_new = cellfun(fun_sqr, cell_core,  'UniformOutput', false);
    core_new = permute(reshape(cell2mat(cell_core_new), 2*(n-1), r0, r1), [2,1,3]);
    cores_new{k} = core_new;
end

% for k = d+1:d+m % last m dimensions for states, compute convolution
%     core = cores_new{k};
%     [r0, ~, r1] = size(core);
%     core_new = core.* eps(Sigma(k-d,k-d)^2 * repmat(freq, r0, 1, r1));
%     cores_new{k} = core_new;
% end

end

