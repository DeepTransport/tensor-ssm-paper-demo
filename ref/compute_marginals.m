function data = compute_marginals(samples, varargin)
% This function prepare the data used for plotting the 1d marginal and 2d
% marginal densities estimated from the samples. 
%
% Note: this function relies on ksdensity function of MATLAB and is 
% expensive to compute. The regula_falsi function for root finding is
% modified from the one in the chebfun package. 
%
% Inputs: 
%
% samples:  dxn matrix, where d is the dimension of the random vector and n
%           is the sample set
%
% Optional parameters:
%
% ngrid:    number of grid points for plotting the marginal densities, the
%           default is 50
%
% qlist:    prescribed quantiles for plotting the contours of 2d marginals,
%           the default is [0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
%
% Output:   data, a struct holding data for plotting
%
% Example:
% 
% load samples.mat
% data = compute_marginals(samples, 'ngrid', 100, 'qlist', [0.9,0.7,0.5,0.3,0.1])
% plot_marginals(data)
%
% Tiangang Cui, August, 2019

defaultNgrid  = 50;
defaultQlist  = [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addRequired (p,'samples',@(x) isnumeric(x));
addParameter(p,'ngrid', defaultNgrid,validScalarPosNum);
addParameter(p,'qlist', defaultQlist,@(x) isnumeric(x));
parse(p,samples,varargin{:});

% samples are supposed to be given as dxn matrix
d = size(samples,1);

nx      = p.Results.ngrid;
qlist   = p.Results.qlist;

data.xi = cell(d,1);
data.fi = cell(d,1);
data.fs = cell(d,d);
data.qants = cell(d,d);

% determine the bounds
minx = zeros(d,1);
maxx = zeros(d,1);
for i = 1:d
    minx(i) = min(samples(i,:));
    maxx(i) = max(samples(i,:));
    tmp     = (maxx(i)-minx(i))/10;
    
    % significant digit
    digi    = 10^floor(log10(tmp));
    l       = floor((minx(i)-tmp*2)/digi)*digi;
    r       = ceil ((maxx(i)+tmp*2)/digi)*digi;
    
    data.xi{i} = linspace(l, r, nx);
end

for i = 1:d
    data.fi{i} = ksdensity(samples(i,:)',data.xi{i},'function','pdf');
    
    % then create two d marginals
    for j = (i+1):d
%         disp([i, j])
        [x1,x2] = meshgrid(data.xi{i},data.xi{j});
        pts     = [x1(:), x2(:)];
        data.fs{i,j}    = ksdensity(samples([i,j],:)', pts, 'BoundaryCorrection', 'Reflection');
        data.qants{i,j} = compute_quantiles(data.fs{i,j}, qlist);
        data.fs{i,j}    = reshape(data.fs{i,j}, size(x1));
    end
end

data.d = d;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qs = compute_quantiles(fs, q_list)

qs   = zeros(size(q_list));

fs   = fs(:);
minf = min(fs);
maxf = max(fs);
z    = sum(fs);
tol  = (maxf - minf)/50;

q0   = minf;
for i = 1:length(q_list)
    qs(i) = regula_falsi(@(x) (sum(fs(fs>x))/z-q_list(i)), tol, q0, maxf);
    q0    = qs(i);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = regula_falsi(func, tol, a, b)


fa = func(a) ;
fb = func(b) ;

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

while ( norm(c-cold, inf) >= tol )
    cold    = c;
    fc  = func(c) ;
    
    I1  = (fc < 0);
    I2  = (fc > 0);
    I3  = ~I1 & ~I2;
    a   = I1.*c + I2.*a + I3.*c;
    b   = I1.*b + I2.*c + I3.*c;
    fa  = I1.*fc + I2.*fa + I3.*fc;
    fb  = I1.*fb + I2.*fc + I3.*fc;
    step    = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
    
    %norm(fc, inf)
    
end

end