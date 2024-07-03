function options = hp_options(varargin)

defaultOptions  = struct([]);
%
defaultMeshSize = 1/40;
defaultXYRatio  = 1;
defaultPolyOrd  = 2;
defaultQuadOrd  = 10;
%
defaultGMRES    = false;
defaultResTol   = 1E-10;
defaultNewtonN  = 1E2;
defaultLinSeaN  = 20;
defaultNLsolver = 'line_search';
expectedNLsolver  = {'line_search', 'trust_region'};
defaultNLdisp   = 1;
%
defaultExpParam = true;
defaultExpThres = 1E-16;
defaultSqParam  = false;
%
defaultObsLoc   = [];
%
defaultBndFuncs = {};
defaultBCTypes  = {};
defaultBCFuncs  = {};
%
defaultQOI      = @NOP;
defaultForce    = @(x1, x2) 0.0;
%
defaultTimeDep  = false;
defaultObsTsta  = inf;
defaultObsTfin  = inf;
defaultObsNtime = 1;
defaultPredT    = inf;
defaultDT       = 0.1;
defaultInit     = @(x1, x2) 0.0;
defaultTlead    = 4;
%

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addOptional (p,'options',   defaultOptions,   @(x) isstruct(x));
addParameter(p,'h',         defaultMeshSize,  validScalarPosNum);
addParameter(p,'xyratio',   defaultXYRatio,   validScalarPosNum);
addParameter(p,'poly_order',defaultPolyOrd,   validScalarPosNum);
addParameter(p,'quad_order',defaultQuadOrd,   validScalarPosNum);
%
addParameter(p,'gmres',     defaultGMRES,     @(x) islogical(x) || isscalar(x));
addParameter(p,'res_tol',   defaultResTol,    validScalarPosNum);
addParameter(p,'num_newton',defaultNewtonN,   validScalarPosNum);
addParameter(p,'num_linsea',defaultLinSeaN,   validScalarPosNum);
addParameter(p,'nl_solver', defaultNLsolver,  @(x) any(validatestring(x,expectedNLsolver)));
addParameter(p,'nl_disp',   defaultNLdisp,    @(x) isnumeric(x) && isscalar(x) && (x >= 0));
%
addParameter(p,'exp_param', defaultExpParam,  @(x) islogical(x) || isscalar(x));
addParameter(p,'exp_thres', defaultExpThres,  validScalarPosNum);
addParameter(p,'sq_param',  defaultSqParam,   @(x) islogical(x) || isscalar(x));
%
addParameter(p,'obs_locs',  defaultObsLoc,    @(x) isnumeric(x));
%
addParameter(p,'bnd_funcs', defaultBndFuncs,  @(x) isa(x, 'cell') );
addParameter(p,'bc_types',  defaultBCTypes,   @(x) isa(x, 'cell') );
addParameter(p,'bc_funcs',  defaultBCFuncs,   @(x) isa(x, 'cell') );
%
addParameter(p,'qoi_func',  defaultQOI,       @(x) isa(x, 'function_handle'));
addParameter(p,'force_func',defaultForce,     @(x) isa(x, 'function_handle'));
%
addParameter(p,'time_dep',  defaultTimeDep,   @(x) islogical(x) || isscalar(x));
addParameter(p,'obs_tstart',defaultObsTsta,   validScalarPosNum);
addParameter(p,'obs_tfinal',defaultObsTfin,   validScalarPosNum);
addParameter(p,'obs_ntime', defaultObsNtime,  validScalarPosNum);
addParameter(p,'pred_t',    defaultPredT,     validScalarPosNum);
addParameter(p,'dt',        defaultDT,        validScalarPosNum);
addParameter(p,'t_lead',    defaultTlead,     validScalarPosNum);
addParameter(p,'init_func', defaultInit,     @(x) isa(x, 'function_handle'));

%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'h','xyratio','poly_order','quad_order','gmres','res_tol',...
        'num_newton','num_linsea','nl_solver','nl_disp','exp_param','exp_thres','sq_param',...
        'obs_locs','bnd_funcs','bc_types','bc_funcs','qoi_func','force_func',...
        'time_dep','obs_tstart','obs_tfinal','obs_ntime','pred_t','dt','init_func','t_lead'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list) || ~isfield(options, list{i})
            options.(list{i}) = tmp.(list{i});
        end
    end
end

if isempty(options.obs_locs)
    error('observation operator is not provided')
end

if options.time_dep
    if isinf(options.obs_tstart) || isinf(options.obs_tfinal) || isinf(options.pred_t)
        error('need to specify time steps for transit problems')
    end
end

end