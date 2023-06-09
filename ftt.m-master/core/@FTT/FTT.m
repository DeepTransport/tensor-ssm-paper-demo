classdef FTT
    % FTT class
    %
    % FTT Properties:
    %   opt         - FTToption
    %   direction   - ALS direction, >0: built from left to right
    %                                <0: built from right to left
    %   oneds       - Data structure containing information for building
    %                 one dimensional polynomial representation.
    %   cores       - Nodal values or coefficent tensors of 1D functions
    %                 the dimension of the current core is organised as:
    %                 previous rank x oned{k}.num_nodes x current rank.
    %   interp_x    - Interpolation coordinates for FTT.
    %   res_x       - Interpolation coordinates for the residual.
    %   res_w       - Interpolation weights for the residual.
    %   n_evals     - Number of function evaluations in TT cross.
    %
    % FTT Methods:
    %   cross       - Run the TT cross, this is called by the constructor.
    %   eval        - Evaluate FTT. The output is horizontally aligned.
    %   eval_block  - Evaluate FTT for either the first or last k variables.
    %                 The output is horizontally aligned.
    %   round       - Round the FTT cores.
    %   int         - Integrate the entire FTT.
    %   int_block   - Integrate a block of FTT cores.
    %   size        - Size of the FTT.
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example 1: (vector function outputs, m = 2):
    %
    % % Step 1: speficy the target function 
    %   func = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
    %   d = 10; % dimensionality of the input
    %
    % % Step 2: setup the Legendre basis polynomial with order 20 in 
    % % domain [0,1]
    %   poly = Legendre(20, [0,1]);
    %
    % % Step 3: use alternating energy enrichment (AMEN), default option
    %   opt = FTToption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, ...
    %           'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
    % % Optional: debug samples
    %   debug_size = 1E4;
    %   debug_x = zeros(d, debug_size);
    %   for k = 1:d
    %       debug_x(k,:) = sample_domain(poly, debug_size);
    %   end
    %
    % % Step 4: build FTT
    %   tt =  FTT(func, d, poly, opt, 'debug_x', debug_x);
    %
    % % Step 5: evaluate the function and its factorisation
    %   exact   = func(debug_x);
    %   appr_tt = eval(tt, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    % % Step 6: round the FTT by truncation local SVD. Here 1E-4 is the 
    % % truncation threshold of each SVD (relative to the largest singular
    % % value)
    %   ttr = round(tt, 1E-4); 
    %   appr_tt = eval(ttr, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %   
    %%%%%%%%%%%%%%%%%
    %
    % Example 2: (use the func and debug samples defined above)
    %
    % % Alternative Lagrange basis function, with order 5 and 4 elements
    %   poly = Lagrangep(5, 4, [0,1], 'ghost_size', 1E-10);
    %
    % % Alternative option use random enrichment
    %   opt_new = update(opt, 'tt_method', 'random');
    %
    % % Build FTT
    %   tt =  FTT(func, d, poly, opt_new, 'debug_x', debug_x);    
    %
    % % Evaluate the function and its factorisation
    %   exact   = func(debug_x);
    %   appr_tt = eval(tt, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also ONED and FTTOPTION
    
    properties
        opt FTToption
        cores
        oneds
        interp_x
        direction
        n_evals
        res_x
        res_w
    end
        
    methods (Static)
        [core, interp_x, res_w, res_x, core_next] = build_basis_amen(oned, ...
            interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr, ...
            dir, int_method, loc_err_tol, max_rank, kick_rank)
        % build tt core by amen

        [core,interp_x,core_next] = build_basis_svd(oned, ...
            interp_xold, core_next, F, ...
            dir, int_method, loc_err_tol, max_rank)
        % build tt core by svd

        [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F)
        % truncate local core

        [f, f_evals] = local_block(sqrt_flag, oned, xleft, xright, func)
        % evaluate function for a coordinate

        interp_x = local_index(oned, direction, interp_xold, ind)
        % index selection

        rerr = local_error(core, f)
        % error of a core
    end
    
    methods
        
        obj = cross(obj, func, d, sample_x, debug_x)
        % cross iterations

        obj = round(obj, thres)
        % round the TT cores
        
        z = int(obj)
        % Integrate the entire TT
        
        ftt = int_block(obj, ind)
        % Integrate a block of TT cores
        
        fx = eval(obj, x, varargin)
        % Evaluate the ftt function. The output is horizontally aligned.
        
        fx = eval_block(obj, x, dir)
        % Evaluate the fTT for either the first or last k variables.
        
        [d,rs,ns] = size(obj)
        % size of ftt

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = FTT(func, d, arg, varargin)
            % Construct tensor train for a function mapping from R^d to R^m.
            %
            %   func - a function (R^d to R^m) that take inputs as a dxn
            %          matrix and returns mxn vector
            %   d    - dimension of the input variable
            %   arg  - either an existing ftt used as the initial guess, or
            %          a set of one dimensional bases for discretising the
            %          function. If only one set of basis is supplied, each
            %          coordinate is discretised using the same basis
            %   opt  - FTT options
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parsing inputs
            defaultOption    = FTToption();
            defaultSampleSet = [];
            defaultDebugSet  = [];
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %
            addRequired(p,'func',@(x) isa(x, 'function_handle'));
            addRequired(p,'d',validScalarPosNum);
            addRequired(p,'arg');
            addOptional(p,'option',defaultOption);
            %
            addParameter(p,'sample_x', defaultSampleSet);
            addParameter(p,'debug_x',  defaultDebugSet);
            %
            p.KeepUnmatched = false;
            parse(p,func,d,arg,varargin{:});
            %
            obj.opt = p.Results.option;
            debug_x = p.Results.debug_x;
            sample_x = p.Results.sample_x;
            %
            if isa(arg, 'FTT')
                obj = arg;
            else
                obj.oneds = cell(d,1);
                if isa(arg, 'cell')
                    for k = 1:d
                        if ~isa(arg{k}, 'oned')
                            error('wrong type of argument')
                        end
                        obj.oneds{k} = arg{k};
                    end
                    obj.cores = [];
                    obj.res_x = [];
                    obj.res_w = [];
                elseif isa(arg, 'oned')
                    for k = 1:d
                        obj.oneds{k} = arg;
                    end
                    obj.cores = [];
                    obj.res_x = [];
                    obj.res_w = [];
                else
                    error('wrong type of argument')
                end
            end
            % build ftt
            obj = cross(obj, func, d, sample_x, debug_x);
        end
        
    end
end