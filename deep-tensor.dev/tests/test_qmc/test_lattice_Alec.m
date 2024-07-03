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
N = 2.^(5:15); % number of QMC samples, changed to powers of 2.
%
Q_qmc = zeros(numel(N), numel(genz)); % RQMC average
Q_mc = zeros(numel(N), numel(genz)); % MC est
%
QMCmse = zeros(numel(N),numel(genz));
MCmse  = zeros(numel(N),numel(genz));

% Initialising random shift
Delta = rand(s, R);

for i = 1:numel(N)
    n = N(i);
    n_mc = R*n;
    %
    Q_r = zeros(R,numel(genz)); % jth randomised QMC rule
    MCest  = zeros(R,numel(genz)); % jth MC block
    MCest2  = zeros(R,numel(genz)); % jth MC block of sum of squares
    for r = 1:R
        Delta_r = Delta(:, r);
        %         % generate all randomly shifted lattice rule samples
        %         X = mod(z*(0:n-1)/n + Delta_r, 1); % each column is a QMC sample
        %         %
        %         Y = rand(s,n); % each column is a MC sample
        %         for j = 1:numel(genz)
        %             Q_j(r,j) = mean(genz{j}(X,w,c));
        %             MCest(r,j) = mean(genz{j}(Y,w,c));
        %             MCest2(r, j) = mean(genz{j}(Y,w,c).^2);
        %         end

        % for loop version
        for k = 0:n-1
            xk = mod(k*z/n + Delta_r, 1); % kth QMC point
            yk = rand(s, 1); % MC sample
            % Computing sum
            for j = 1:numel(genz)
                Q_r(r,j) = Q_r(r, j) + genz{j}(xk,w,c);
                MCest(r,j) = MCest(r, j) + genz{j}(yk,w,c);
                MCest2(r, j) = MCest2(r, j) + genz{j}(yk,w,c).^2;
            end
        end
        % Averaging
        Q_r(r, :) = Q_r(r, :)/n;
        MCest(r, :) = MCest(r, :)/n;
        MCest2(r, :) = MCest2(r, :)/n;
    end
    % Averaging over all randomisations/MC blocks
    Q_qmc(i, :) = mean(Q_r, 1);
    Q_mc(i, :) = mean(MCest, 1);
    % Estimating MSE by sample variance
    QMCmse(i,:) = sum((Q_r - Q_qmc(i, :)).^2, 1)/R/(R - 1);
    MCmse(i,:) = sum((MCest - Q_mc(i, :)).^2, 1)/R/(R - 1);
    %MCmse(i,:) = (mean(MCest2, 1) - Q_mc(i, :).^2)/(n_mc - 1);
end

%%
for j = 1:numel(genz)
    fig = figure(j);
    loglog(N, sqrt(MCmse(:,j)), '-r*')
    hold on
    plot(N, sqrt(QMCmse(:,j)), '-bo')
    plot(N, 1./sqrt(N)/100, '--r')
    plot(N, 1./N/100, '--b')
    hold off
    title(['genz ' num2str(j)])
    xlabel('N')
    ylabel('rMSE')
    h = legend('MC', 'QMC', '$\sim 1/\sqrt{N}$', '$\sim 1/N$');
    set(h, 'Interpreter', 'latex', 'fontsize', 16)
end

