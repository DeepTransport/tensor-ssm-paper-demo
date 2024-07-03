classdef SparseFun < ApproxFun
    
    properties (Constant)
        defaultVar = InputData()
        defaultOpt = SparseOption()
        defaultData = SparseData()
        defaultPoly = Legendre(30)
        defaultDomain = BoundedDomain([-1,1])
    end
    
    methods (Static)
        [x,err,Q,R] = ls_solve(A, y, w)
        [x,err,Q,R] = ls_solve_update(Q,R,A,y,w)
    end
    
    properties
        importance % Christoffel function
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function n = cardinal(obj)
            n = cardinal(obj.data.I);
        end
        
        function flag = isweighted(obj)
            switch obj.opt.weight_rule
                case {'QMC'}
                    flag = false;
                case {'Christoffel'}
                    flag = true;
                otherwise
                    flag = false;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function x = new_QMC_points(obj, log2N, shift)
            u = new_points(obj.importance, log2N, shift);
            %u = u(:,2:end);
            x = measure_reference_inverse_cdf(obj.base, u);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x = add_QMC_points(obj, log2N, shift)
            u = add_points(obj.importance, log2N, shift);
            x = measure_reference_inverse_cdf(obj.base, u);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function x = random(obj, I, sample_factor)
            if isweighted(obj)
                x = random(obj.importance, I, sample_factor);
            else
                x = sample_measure_reference(obj.base, cardinal(I)*ceil(sample_factor));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function w = eval_weight(obj, I, x)
            if isweighted(obj)
                w = eval_weight(obj.importance, I, x);
                w = w(:);
            else
                w = [];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function basis_at_z = eval_basis(obj, I, z)
            % basis_at_z: num(x) x cardinal
            if isa(I, 'MultiIndices')
                I = I.array;
            end
            basis_at_z = ones(size(z,2), size(I,1));
            for k = 1:size(I,2)
                ind = I(:,k);
                bk = eval_basis(obj.base.oneds{k}, reshape(z(k,:),[],1));
                basis_at_z = basis_at_z.*bk(:,ind+1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_reference(obj, z)
            basis_at_z = eval_basis(obj, obj.data.I, z);
            f = reshape(basis_at_z*obj.data.coeff,[],size(z,2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [g,f] = grad_reference(obj, z)
            basis_at_z = eval_basis(obj, obj.data.I, z);
            g = zeros(ndims(obj), size(z,2));
            for k = 1:ndims(obj)
                ind = obj.data.I.array(:,k);
                bk  = eval_basis(obj.base.oneds{k}, reshape(z(k,:),[],1));
                dbk = eval_basis_deri(obj.base.oneds{k}, reshape(z(k,:),[],1));
                dbk = dbk./bk;
                dbasis_at_z = basis_at_z.*dbk(:,ind+1);
                g(k,:) = reshape(dbasis_at_z*obj.data.coeff,[],size(z,2));
            end
            f = reshape(basis_at_z*obj.data.coeff,[],size(z,2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = int_reference(obj)
            z = obj.data.coeff(1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = SparseFun(func, arg, varargin)
            obj@ApproxFun(func, arg, varargin{:});
            %
            switch obj.opt.weight_rule
                case {'QMC'}
                    obj.importance = Lattice(ndims(obj.base));
                    obj = least_square_qmc(obj, func);
                case {'Christoffel'}
                    obj.importance = Christoffel(obj.base);
                    obj = least_square_mc(obj, func);
                otherwise
                    obj.importance = [];
                    obj = least_square_mc(obj, func);
            end
        end
    end
    
end