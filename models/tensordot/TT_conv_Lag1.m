function yval = TT_conv_Lag1(ftt, Sigma_update, x)
% Input is a tensor train ,sigma_update function, and x values
% Output is the function evaluation of the convolved function at x
% The ftt is established wrt Lagrange1 function

% m = size(Sigma, 1);

Sigma = Sigma_update(x); % Sigma is of size m * nx
dm = size(x, 1);
nx = size(x, 2);


% ind = (any(x <  -1)) | ( any(x > 1));

% order = (size(ftt.data.cores{1}, 2) - 1)/2 ;
% freq2 = [0, 1:order, 1:order].^2;
% 
% % modifier(dm-m+1:dm, :,:) = exp(-0.5*pi^2*repmat(diag(Sigma), 1, nx).* reshape(freq2, 1,1,[]));
% modifier = exp(-0.5*pi^2* Sigma.* reshape(freq2, 1,1,[]));
frl  = ones(nx,1);
% fradon = ones(1,nx);
x_ref = domain2reference(ftt.base.oned_domains{1}, x);
for j = 1:min(dm,ndims(ftt.base))
    %
    % pk  = abs(frl*obj.ys{j})';
    % %
    % tmp = eval_pdf_radon(obj.oned_cdfs{j}, pk+obj.tau, x_ref(j,:));
    % fradon = fradon.*reshape(tmp,size(fradon));

    nj  = size(ftt.data.cores{j}, 2);
    rjm = size(ftt.data.cores{j}, 1);
    %
    if j < ndims(ftt.base) || (size(ftt.data.cores{j}, 3) > 1 && size(ftt.data.cores{j}, 4) == 1)
        tmp = reshape(permute(ftt.data.cores{j}, [2,1,3]), nj, []);
    else
        % collapse the third dimension = 1
        tmp = reshape(permute(reshape(ftt.data.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
    end
    domainx = reference2domain(ftt.base.oned_domains{1}, ftt.base.oneds{1}.nodes);
    T   = reshape(permute(reshape(eval_radon_conv(ftt.base.oneds{1}, tmp, x(j,:), x_ref(j, :), Sigma(j, :), domainx), nx, rjm, []), [2,1,3]), rjm*nx, []);
    % how to speed up this part?
    jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rjm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rjm*nx);
    %
    frl  = B*T;
end
yval = frl';


    function f = eval_radon_conv(obj, coeff, x, x_ref, Sigma, domainx)
        if sum(Sigma == 0) == sum(Sigma)
            f = eval_basis(obj, x_ref')*coeff;
        else
            dx = domainx(2) - domainx(1);
            Sigma = sqrt(Sigma);
            xminusover = (x - domainx)./(Sigma.*sqrt(2)); % l * nx matrix
            ell = size(coeff, 1);
    
            c1 = reshape(diff(coeff)/sqrt(2*pi)/dx, ell-1, 1, []);
            c2 = reshape(c1*sqrt(pi/2), ell-1, 1, []);
            c3 = reshape(0.5*(domainx(2:end).*coeff(1:end-1, :) - domainx(1:end-1).*coeff(2:end, :))/dx, ell-1, 1, []);
    
            temp = erf(xminusover(1:end-1, :)) - erf(xminusover(2:end, :));
            M1 = c1 .* (exp(-xminusover(1:end-1, :).^2) - exp(-xminusover(2:end, :).^2)).*Sigma;
            M2 = c2 .* temp .* x;
            M3 = c3 .* temp;
    
            f = reshape(sum(M1+M2+M3), size(x, 2), []);
            f(Sigma == 0, :) = eval_basis(obj, x_ref(Sigma == 0)')*coeff;
        end
        ind = (x < domainx(1)) | (x > domainx(end));
        f(ind, :) = 0;
    end

end

