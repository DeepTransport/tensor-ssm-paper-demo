% % test precond
% m = 4;
% d = 3;
% A = randn(m+d,m+d);
% A = A*A';
% B = randn(m,m);
% B = B*B';
% B2 = [zeros(d), zeros(d,m);zeros(m,d), B];
% [C, D] = precond(A, B);
% 
% C*A*C'
% C(d+1:d+m,d+1:d+m)*B*C(d+1:d+m,d+1:d+m)'
% C*B2*C'
% 
% % precond pass
% 
% % test cellfun
% n = 2;
% A = rand(n,6);
% B = rand(n,6);
% A + B
% A = mat2cell(A, n, ones(6,1));
% B = mat2cell(B, n, ones(6,1));
% 
% C = cellfun(@plus, A, B,  'UniformOutput', false);
% cell2mat(C)
% % cellfun pass
%% test eval_radon_conv
% dom = BoundedDomain([-4,4]);
% % poly = bases{1};
% poly = ApproxBases(Lagrange1(1000), dom, dd);
% domainx = linspace(-4,4,1001)';
% fun = @(x) normpdf(x, 0, sqrt(0.9));
% coeff = fun(domainx);
% Sigma = 0.1;
% fun_conv = @(x) normpdf(x, 0, sqrt(0.9+Sigma));
% 
% % x = randn(1, 20);
% % x = 0.1;
% x = domainx';
% y_true = fun_conv(x);
% c = 20;
% y_app = eval_radon_conv(poly.oneds{1}, c*coeff, x, x/4, Sigma, domainx)'/c;
% 
% y_true./y_app
% 
% hold  on
% plot(y_true)
% plot(y_app)
% hold off
%%
% test TT_sqrconv
dd = 3 ;
dom = BoundedDomain([-4,4]);
bases{1} = ApproxBases(Fourier(20), dom, dd);
bases{2} = ApproxBases(Fourier(40), dom, 2);
bases{3} = ApproxBases(Lagrange1(40), dom, 2);
bases{4} = ApproxBases(Lagrangep(5,8), dom, 2);
% poly = bases{1};
poly = ApproxBases(Lagrange1(100), dom, dd);

opts{1} = TTOption('tt_method', 'amen', ...
    'als_tol', 1E-2, 'local_tol', 1E-2, 'max_rank', 30, 'max_als', 6, 'init_rank', 20, 'kick_rank', 5);
opts{2} = TTOption('tt_method', 'random', ...
    'als_tol', 1E-8, 'local_tol', 1E-20, 'max_rank', 19, 'max_als', 6);
opt = opts{1};
sigma = diag([0.1,0.2,0.3]);
% sigma = zeros(3);
% M = poly.oneds{1}.basis2node;

x = repmat(linspace(-4, 4, 100), 3,1);

% fun = @(x) (3+sin(pi * x(1,:))).^2.*(2+cos(pi* x(2,:))).^2;
fun = @(x) mvnpdf(x', 0, eye(dd))';
fun_conv = @(x) mvnpdf(x', 0, eye(dd)+ sigma)';
% C = randn(dd);
% fun = @(x) 0.5*sum((C*x).^2,1);

sirt = TTIRT(@(x) -log(fun(x)), poly, opt);
eval_pdf(sirt, [1;1;1])

fun(x)./eval_pdf(sirt, x)
figure
hold on
plot(fun(x))
plot(eval_pdf(sirt, x))
hold off

fun1 = @(x) mvnpdf(x', 0, .5*eye(dd))';
sirt1 = TTSIRT(@(x) -log(fun1(x)), poly, opt);
eval_pdf(sirt1, [1;1;1])
% fun(x)./eval_pdf(sirt1, x)
% figure
% hold on
% plot(fun(x))
% plot(eval_pdf(sirt1, x))
% hold off

% eval_pdf(sirt1,x)
% sirt2 = TTIRT(@(x) -0.5*log(fun(x)), poly, opt);

ftt = sirt.approx;
% sirt  = TTSIRT(fun, poly, opt); % , 'samples', C\rand(dd, 100)

fun_updatecov = @(x) diag(sigma).*ones(1, size(x, 2));
% domainx = reference2domain(ftt.base.oned_domains{1}, ftt.base.oneds{1}.nodes);


% x = [0;0;0];
% eval_pdf(sirt, x)
y_true = fun_conv(x);
y_true = y_true/sum(y_true);
y_ftt = TT_conv_Lag1(ftt, fun_updatecov, x);
y_ftt = y_ftt/sum(y_ftt);
y_true./abs(y_ftt)
% eval_pdf(sirt, rand(2,5))

figure
hold on
plot(y_true)
plot(y_ftt)
hold off



%% test ttdot
r = [1,3,5,6,1];
s = [1,6,7,8,1];
n = [2,5,6,7];
cores1 = cell(4);
cores2 = cell(4);

for k = 1:4
    cores1{k} = randn(r(k), n(k), r(k+1));
    cores2{k} = randn(s(k), n(k), s(k+1));
end

a1 = zeros(n);
a2 = zeros(n);

for k1 = 1:n(1)
    for k2 = 1:n(2)
        for k3 = 1:n(3)
            for k4 = 1:n(4)
                a1(k1,k2,k3,k4) = reshape(cores1{1}(1,k1,:), 1, [])*reshape(cores1{2}(:,k2,:), r(2), [])...
                    *reshape(cores1{3}(:,k3,:), r(3), [])*reshape(cores1{4}(:,k4,1), r(4), []);
                a2(k1,k2,k3,k4) = reshape(cores2{1}(1,k1,:), 1, [])*reshape(cores2{2}(:,k2,:), s(2), [])...
                    *reshape(cores2{3}(:,k3,:), s(3), [])*reshape(cores2{4}(:,k4,1), s(4), []);
            end
        end
    end
end

c_full = sum(a1.*a2, "all");
c_tt = ttdot(cores1, cores2);


% ttdot pass


%%
% test TT_eval
dd = 8;
C = spdiags([0.5*ones(dd,1), ones(dd,1), 0.5*ones(dd,1)], [-1,0,1],dd,dd);
C = (C+C')/2;
D = diag([0,0, 2*rand(1,dd-2)]);
% diag(2*rand(1,dd));
% D = diag([1,2,3:dd]); 
% D = zeros(dd);
% D = 4e-9*eye(dd);


dom = BoundedDomain([-4,4]);
bases{1} = ApproxBases(Fourier(5), dom, dd);
bases{2} = ApproxBases(Fourier(20), dom, dd);
bases{3} = ApproxBases(Lagrange1(40), dom, 2);
bases{4} = ApproxBases(Lagrangep(5,8), dom, 2);
poly = bases{2};
poly2 = ApproxBases(Fourier(40), dom, dd);


opts{1} = TTOption('tt_method', 'amen', ...
    'als_tol', 1E-2, 'local_tol', 1E-2, 'max_rank', 30, 'max_als', 8, 'init_rank', 10, 'kick_rank', 5);
opts{2} = TTOption('tt_method', 'random', ...
    'als_tol', 1E-8, 'local_tol', 1E-20, 'max_rank', 19, 'max_als', 6);
opt = opts{1};



fun_origin = @(x) mvnpdf(x', zeros(1,dd), C)';
fun_conv = @(x) mvnpdf(x', zeros(1,dd), C+D)';

% rng(0)
sirt = TTSIRT(@(x) -log(fun_origin(x)), poly, opt);
for k = 1:dd
    sirt.approx.data.cores{k}(:, end, :) = 0;
end


z = chol(C)*randn(dd, 10);
y1 = fun_origin(z);
fun_ftt = @(x) eval_pdf(sirt, x);
y2 = fun_ftt(z);
disp(y1./y2)
% sirt  = TTSIRT(fun, poly, opt); % , 'samples', C\rand(dd, 100)

% [tt_new, ploy_new]= TT_sqrconv(sirt, sigma);
% eval_pdf(sirt, rand(2,5))
TT_new = TT_expand(sirt);

fun3 = @(x) TT_eval(TT_new, poly2, @(x) diag(D).*ones(1, size(x, 2)), x);

% z = randn(dd, 10);
y3 = fun_conv(z);
y4 = fun3(z);
disp(y3./y4)
rng('shuffle')

sirt2 = sirt;
TT_sqr_cores = cell(dd, 1);
for k = 1:dd
   [r0, l, r1] = size(TT_new{k});
   corenew = zeros(r0,l,r1);
   for k1 =  1:r0
      for k2 = 1:r1
          corenew(k1, :, k2) = poly2.oneds{1}.node2basis*TT_new{k}(k1, :, k2)';
      end
   end
   TT_sqr_cores{k} = corenew;
end

sirt2.approx.data.cores = TT_sqr_cores;
sirt2.approx.base.oneds = poly2.oneds;
sirt2 = marginalise(sirt2, 1);
fun_sqr = @(x) eval_pdf(sirt2, x);
y_sqr = fun_sqr(z);
disp(y1./y_sqr);

function f = eval_radon_conv(obj, coeff, x, x_ref, Sigma, domainx)
        dx = domainx(2) - domainx(1);
        Sigma = sqrt(Sigma);
        xminusover = (x - domainx)./(Sigma.*sqrt(2)); % l * nx matrix

        c1 = diff(coeff)/sqrt(2*pi)/dx;
        c2 = c1*sqrt(pi/2);
        c3 = 0.5*(domainx(2:end).*coeff(1:end-1) - domainx(1:end-1).*coeff(2:end))/dx;

        temp = erf(xminusover(1:end-1, :)) - erf(xminusover(2:end, :));
        M1 = c1 .* (exp(-xminusover(1:end-1, :).^2) - exp(-xminusover(2:end, :).^2)).*Sigma;
        M2 = c2 .* temp .* x;
        M3 = c3 .* temp;

        f = sum(M1+M2+M3)';
        f(Sigma == 0) = eval_basis(obj, x_ref(Sigma == 0)')*coeff;
        ind = (x < domainx(1)) | (x > domainx(end));
        f(ind,:) = 0;
end
