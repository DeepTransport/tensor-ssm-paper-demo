function output = ftt_sqrconv(ftt, poly1, poly2, d, sigma2)
% compute the square of ftt and then convolution with independent Gaussian
% sigma2 is a column vector storing the variance of eahc dim of Gaussian
% ftt is the FTT for the sqrt of the target function
% d is the dimension for theta, first d dimensions no convolution
% arguments
%     ftt1 FTT
%     ftt2 FTT
% end
output = ftt;


% check important properties of ftt1 and ftt2 are the same
% if ftt1.direction ~= ftt2.direction, error('the ftt directions not match'); end
% if ftt1.direction ~= ftt2.direction, ftt2 = marginalise(ftt2, ftt1.direction); end
% ftt2 = marginalise(ftt2, ftt1.direction);


output.interp_x = [];
output.l2_err = [];
output.n_evals = [];
output.res_x = [];
output.res_w = [];

cores_ftt = ftt.cores;
cores_out = cores_ftt;


c = ftt1.oned_domains{1}.bound(2);
freq = ((1:poly2.m - 1)/c)';
FG = exp(-0.5*pi^2*freq.^2*sigma2');
FG_all = [ones(1, size(FG, 2)); repmat(FG,2,1)];
FG_all = [ones(size(FG_all, 1), d), FG_all];


% build matrix FG first


for k = 1:size(cores_ftt, 1) % for each x_k
    core = cores_ftt{k};
    r1 = size(core, 1);
    r2 = size(core, 3);
    core_out = zeros(r1^2, size(core, 2), r2^2);
    CL = cell(r1^2, r2^2);
    ind1r = kron(ones(r1,r2^2), 1:r1);
    ind1c = kron(ones(r1^2, r2), 1:r2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% here!!!!!!!!!!!!!!!!!!!!

    for j = 1:size(core1, 2) % for each dimension of basis function, from 1 to n_k
        temp = kron(reshape(core1_exp(:,j,:), size(core1_exp(:,j,:),1), size(core1_exp(:,j,:), 3)),...
            reshape(core2_exp(:,j,:), size(core2_exp(:,j,:),1), size(core2_exp(:,j,:), 3)));
        core_out_exp(:,j,:) = temp;
    end
    cores_out{d + k} = real(change_basis(core_out_exp, M, -1));
end
output.cores = cores_out;
disp("rounding for the convolution")
output = round(output, 1e-10);

end