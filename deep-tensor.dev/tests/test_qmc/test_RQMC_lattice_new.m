%%
s = 10; % dimension

% define genz functions
w = rand(1,s)*0.6 + 0.2;
c = rand(1,s)*0.6 + 0.2;
%
genz{1} = @(x, w, c) cos(2*pi*w(1) + c(:)'*x);
genz{2} = @(x, w, c) exp( - c(:)'*(x-w(:)).^2 );
genz{3} = @(x, w, c) exp( - c(:)'*abs(x-w(:)) );

%
R = 32; % number of randomised QMC sample sets

% Getting lattice rule
fid = fopen('lattice-39101-1024-1048576.3600');
[z, ~] = fscanf(fid, '%d %d', [2, inf]);
z = z(2, 1:s)';

%%
log2N = 5:18; % number of QMC samples, changed to powers of 2.
%
Q_qmc = zeros(numel(genz),numel(log2N)); % RQMC estimates
Q_mc  = zeros(numel(genz),numel(log2N)); % MC estimates
%
QMCmse = zeros(numel(genz),numel(log2N));
MCmse  = zeros(numel(genz),numel(log2N));

% Initialising random shift
Delta = rand(s, R);

%
ll = Lattice(s);

sum_qmc = zeros(numel(genz),R); % RQMC outputs
sum_mc  = zeros(numel(genz),R); % MC outputs
for i = 1:numel(log2N)
    n = 2^log2N(i);
    for r = 1:R
        if i == 1
            QMC_shifted = new_points(ll,log2N(i),Delta(:,r));
            N_new = 2^log2N(i);
        else
            QMC_shifted = add_points(ll,log2N(i-1),Delta(:,r));
            N_new = 2^log2N(i-1);
        end
        MC_new = rand(s,N_new);
        for j = 1:numel(genz)
            sum_qmc(j,r) = sum_qmc(j,r) + sum(genz{j}(QMC_shifted,w,c));
            sum_mc(j,r)  = sum_mc(j,r) + sum(genz{j}(MC_new,w,c));
        end
    end
    % Averaging
    E_qmc = sum_qmc/n;
    E_mc  = sum_mc/n;
    % Averaging estimator
    Q_qmc(:,i) = mean(E_qmc,2);
    Q_mc(:,i) = mean(E_mc,2);
    % Estimating MSE by sample variance
    QMCmse(:,i) = sum((E_qmc - Q_qmc(:,i)).^2, 2)/R/(R - 1);
    MCmse(:,i)  = sum((E_mc - Q_mc (:,i)).^2, 2)/R/(R - 1);
end

%%
for j = 1:numel(genz)
    fig = figure(j);
    loglog(2.^log2N, sqrt(MCmse(j,:)), '-r*')
    hold on
    plot(2.^log2N, sqrt(QMCmse(j,:)), '-bo')
    plot(2.^log2N, 1./sqrt(2.^log2N)/100, '--r')
    plot(2.^log2N, 1./(2.^log2N)/100, '--b')
    hold off
    title(['genz ' num2str(j)])
    xlabel('N')
    ylabel('rMSE')
    h = legend('MC', 'QMC', '$\sim 1/\sqrt{N}$', '$\sim 1/N$');
    set(h, 'Interpreter', 'latex', 'fontsize', 16)
end

