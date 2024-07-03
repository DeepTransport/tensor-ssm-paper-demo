function options = inverse_options(varargin)
    
defaultOptions  = struct([]);
%
%
defaultStd      = 1;
defaultS2n      = 10;
defaultJacobian = false;
defaultLike     = 'normal';
expectedLike    = {'normal', 'poisson'};
%
defaultCovType  = 'MRF';
expectedCovType = {'Conv', 'GP', 'MRF'};
defaultMean     = 0;
defaultScale    = 0.25*eye(2);
defaultPower    = 1;
defaultSigma    = 1;
defaultK        = 1E-15;
defaultCond     = [1, 1, 0];
defaultCenters  = [];
defaultRadii    = [];
defaultWeights  = [];
defaultKernel   = @(d) exp( - 0.5*d.^2);
defaultTestType  = 'Prior';
expectedTestType = {'CF', 'Inclusion', 'Prior'};
defaultTestBase  = 2;
defaultTestRange = 2;

%
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
%
addOptional (p,'std',       defaultStd,     validScalarPosNum);
addParameter(p,'s2n',       defaultS2n,     validScalarPosNum);
addOptional (p,'jacobian',  defaultJacobian,@(x) islogical(x) || isscalar(x));
addParameter(p,'like',      defaultLike,    @(x) any(validatestring(x,expectedLike)));
%
addOptional (p,'options',   defaultOptions, @(x) isstruct(x));
%
addParameter(p,'cov_type',  defaultCovType, @(x) any(validatestring(x,expectedCovType)));
addOptional (p,'mean',      defaultMean,    @(x) isnumeric(x));
addOptional (p,'scale',     defaultScale,   @(x) isnumeric(x));
addOptional (p,'power',     defaultPower,   @(x) isnumeric(x));
addOptional (p,'sigma',     defaultSigma,   validScalarPosNum);
addOptional (p,'gamma',     defaultK,       validScalarPosNum);
addOptional (p,'k',         defaultK,       validScalarPosNum);
addOptional (p,'cond',      defaultCond,    @(x) isnumeric(x));
addOptional (p,'centers',   defaultCenters, @(x) isnumeric(x));
addOptional (p,'radii',     defaultRadii,   @(x) isnumeric(x));
addOptional (p,'weights',   defaultWeights, @(x) isnumeric(x));
addParameter(p,'kernel_func', defaultKernel,@(x) isa(x, 'function_handle'));
%
addParameter(p,'test_type', defaultTestType, @(x) any(validatestring(x,expectedTestType)));
addOptional (p,'test_base', defaultTestBase);
addOptional (p,'test_range',defaultTestRange);
%
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'std', 's2n', 'jacobian', 'like',...
        'cov_type', 'mean', 'scale', 'power', 'sigma', 'gamma', 'k', 'cond', ...
        'centers', 'radii', 'weights', 'kernel_func', 'test_type', 'test_base', 'test_range'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list)
            options.(list{i}) = tmp.(list{i});
        end
    end
end

if strcmp(options.cov_type, 'Conv')
    nc = size(options.centers,1);
    if nc < 1
        error('needs more than 1 kernel centers')
    end
    if nc ~= length(options.radii) || nc ~= length(options.weights)
        error('number of centers needs to math number of radii and number of weights')
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

