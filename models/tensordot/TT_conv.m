function yval = TT_conv(ftt, Sigma_update, x)
% Input is a tensor train ,sigma_update function, and x values
% Output is the function evaluation of the convolved function at x

% m = size(Sigma, 1);

Sigma = Sigma_update(x);
dm = size(x, 1);
nx = size(x, 2);


% ind = (any(x <  -1)) | ( any(x > 1));

order = (size(ftt.cores{1}, 2) - 1)/2 ;
freq2 = [0, 1:order, 1:order].^2;

% modifier(dm-m+1:dm, :,:) = exp(-0.5*pi^2*repmat(diag(Sigma), 1, nx).* reshape(freq2, 1,1,[]));
modifier = exp(-0.5*pi^2* Sigma.* reshape(freq2, 1,1,[]));
fx  = ones(nx,1);
% all the intermediate dimensions, except the last dimension

for j = 1:min(dm,ndims(ftt))
    nj  = size(ftt.cores{j}, 2);
    rjm = size(ftt.cores{j}, 1);
    %
    if j < ndims(ftt) || (size(ftt.cores{j}, 3) > 1 && size(ftt.cores{j}, 4) == 1)
        tmp = reshape(permute(ftt.cores{j}, [2,1,3]), nj, []);
    else
        % collapse the third dimension = 1
        tmp = reshape(permute(reshape(ftt.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
    end
    T   = reshape(permute(reshape(eval_radon_conv(ftt.oneds{j},tmp,x(j,:), reshape(modifier(j,:,:), nx, []) ), nx, rjm, []), [2,1,3]), rjm*nx, []);
    % how to speed up this part?
    jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rjm);
    B   = sparse(ii(:), jj(:), fx(:), nx, rjm*nx);
    %
    fx  = B*T;
end

yval = fx';
% yval(1, ind) = 0;

    function f = eval_radon_conv(obj, coeff, x, modifier)
        % Evaluates the approximated function for a given vector of x
        % points. f = EVAL(poly, A, x)
        %
        %   A - Either the nodal values (piecewise class) or the
        %       coefficients (spectral class), dim(basis) x num(x)
        %   x - A vector of x points.
        %   f - A column vector of outputs num(x) x 1

        b = eval_basis(obj, x(:)).*modifier;
        f = b*coeff;
        %
        ind = (x < obj.domain(1)) | (x > obj.domain(2));
        f(ind,:) = 0;
    end

end

