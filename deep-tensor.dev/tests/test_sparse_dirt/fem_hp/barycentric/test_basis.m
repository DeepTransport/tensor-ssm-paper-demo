function test_basis(num_dim, poly_order, quad_order, m, n)
%test_basis
%
% test the basis function 
%
% Tiangang Cui, 29/Oct/2016

[local_elem, local_elem_bnd] = local_elem_simplex(num_dim, poly_order, quad_order);

disp('quad elem')
test_quad(local_elem, m);

disp('quad elem bnd')
test_quad(local_elem_bnd, m);

disp('grad elem')
test_grad(local_elem, n);

disp('grad elem bnd')
test_grad(local_elem_bnd, n);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_quad(elem, n)

syms a b c


ai = zeros(n,1);
qi = zeros(n,1);
ii = zeros(n,1);

if elem.num_dim == 1
    c = rand(n, 1);
    r = rand(n, 1);
    p = rand(n, 1);
    p = ceil(p*elem.quad_order/2);
    for i = 1:n
        f(a) = c(i) + r(i)*a^p(i);
        ai(i) = eval(int(f, a, 0, 1));
        qf = elem.f_quad_pts*eval(f(elem.lambda_pts));
        ii(i) = sum(elem.quad_weights.*qf);
        qi(i) = sum(elem.quad_weights.*eval(f(elem.quad_lambdas)));
    end
end

if elem.num_dim == 2
    c = rand(n, 1);
    r = rand(n, 3); 
    p = rand(n, 4);
    p(:,1:2) = floor(p(:,1:2)*elem.quad_order/2);
    p(:,3:4) = floor(p(:,3:4)*elem.quad_order/4);
    for i = 1:n
        f(a, b) = c(i) + r(i,1)*a^p(i,1) + r(i,2)*b^p(i,2) + r(i,3)*a^p(i,3)*b^p(i,4);
        ai(i) = int( int(f, b, 0, 1-a), a, 0, 1 );
        qf = elem.f_quad_pts*eval(f(elem.lambda_pts(:,1), elem.lambda_pts(:,2)));
        ii(i) = sum(elem.quad_weights.*qf);
        qi(i) = sum(elem.quad_weights.*eval(f(elem.quad_lambdas(:,1), elem.quad_lambdas(:,2))));
    end
end

%[ai, qi] 

disp(mean(abs(ai(:) - ii(:)))/mean(abs(ai(:))))
    
disp(mean(abs(ai(:) - qi(:)))/mean(abs(ai(:))))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_grad(elem, n)

% generate a set of points
tol = 1E-6;
r = rand(n, elem.num_dim+1);
lambdas = r(:, 1:end-1)./repmat(sum(r, 2), 1, elem.num_dim);

for d = 1:elem.num_dim
    % perturb 
    lambdas_p = lambdas;
    lambdas_m = lambdas;
    lambdas_p(:,d) = lambdas_p(:,d) + tol;
    lambdas_m(:,d) = lambdas_m(:,d) - tol;
    
    f_p = eval_barycentric_bases(elem.alpha, elem.A, lambdas_p);
    f_m = eval_barycentric_bases(elem.alpha, elem.A, lambdas_m);
    
    gd  = (f_p - f_m)/(2*tol);
    g   = eval_barycentric_bases(elem.grad_alpha{d}, elem.grad_A{d}, lambdas);
    
    disp(['grad ' num2str(d)])
    disp(mean(abs(gd(:) - g(:)))/mean(abs(g(:))))
%     figure
%     subplot(2,1,1)
%     plot(abs(gd(:) - g(:)))
%     title(['grad ' num2str(d)])
%     subplot(2,1,2)
%     plot(gd(:), g(:))
end

end