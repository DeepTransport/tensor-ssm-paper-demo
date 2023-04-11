classdef oned
    % oned class - Superclass for all one dimensional basis
    %
    % oned Properties:
    %   domain          - The domain of discretisation.
    %   order           - Order of the discretisation basis.
    %   num_nodes       - Number of collocation points.
    %   nodes           - A column vector of ollocation points.
    %
    % oned Methods:
    %   eval            - Evaluates the approximated function for a given  
    %                     vector of x points.
    %   eval_deri       - Evaluates the derivative of the approximated 
    %                     function for a given vector of x.
    %   eval_basis      - Evaluates one dimensional basis on a given vector
    %                     of x points. num(x) x dim(basis)
    %   eval_basis_deri - Evaluates the derivative of the one dimensional
    %                     basis for a given vector of x. num(x) x num(basis)
    %   sample_domain   - Generates random samples from the one dimensional
    %                     approximation domain.
    %   mass_r          - Evaluates the matrix vector (or matrix) product 
    %                     with the upper Cholesky factor of the mass matrix. 
    %                     Used in the squared Rosenblatt transport. 
    %   integral        - Evaluates the oned dimensional integral given 
    %                     either evaluation at collocation points (piecewise 
    %                     basis) or the coefficents of spectral basis. 
    %   point_selection - Select interpolation points for constructing FTT.
    %
    %%%%%%%%%%%%%%%%%
    %
    % Two subclasses are implemented. The piecewise class implements bases
    % using nodal values, which contains piecewise linear basis (Lagrange1)
    % and piecewise high order basis (Lagrangep). The spectral class has
    % spectral bases, including spectral polynomials, under the recurr
    % subclass (Legendre, Jabobi, Hermite, and Laguerre), Chebyshev 1st and 
    % 2nd polynomials subclasses, and Fourier subclass.  
    %
    % See also PIECEWISE, SPECTRAL, and ONEDCDF
    
    properties
        domain(1,2)
        order
        num_nodes
        nodes(:,1)
        mass_R(:,:)
        int_W(1,:)
    end
    
    methods (Abstract)
        eval_basis(obj)
        eval_basis_deri(obj)
        %
        sample_domain(obj)
        point_selection(obj)   
    end
    
    methods
        function f = eval(obj, coeff, x)
            % Evaluates the approximated function for a given vector of x 
            % points. f = EVAL(poly, A, x)
            %   
            %   A - Either the nodal values (piecewise class) or the
            %       coefficients (spectral class), dim(basis) x num(x)
            %   x - A vector of x points. 
            %   f - A column vector of outputs num(x) x 1
            
            
            b = eval_basis(obj, x(:));
            f = b*coeff;
            ind = (x < obj.domain(1)) | (x > obj.domain(2));
            f(ind,:) = 0;
            %{
            f = zeros(length(x), size(coeff,2));
            mid = (x >= obj.domain(1)) & (x <= obj.domain(2));
            b = eval_basis(obj, x(mid));
            f(mid,:) = b*coeff;
            %}
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_deri(obj, coeff, x)
            % Evaluates the derivative of the approximated function for a 
            % given vector of x points. f = EVAL_DERI(poly, A, x)
            %   
            %   A - Either the nodal values (piecewise class) or the
            %       coefficients (spectral class), dim(basis) x num(x)
            %   x - A vector of x points. 
            %   f - A column vector of outputs num(x) x 1

            f   = zeros(length(x), size(coeff,2));
            mid = (x >= obj.domain(1)) & (x <= obj.domain(2));
            %
            if sum(mid) > 0
                b = eval_basis_deri(obj, reshape(x(mid),[],1));
                f(mid,:) = b*coeff;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function R = mass_r(obj, interp_w)
            % Evaluates the matrix vector (or matrix) product with the 
            % upper Cholesky factor of the mass matrix. 
            %   R = MASS_R(poly, interp_w)
            %
            %   interp_w    - Either the nodal values (piecewise) or 
            %                 corefficents (spectral)
            
            R = obj.mass_R*interp_w;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = integral(obj, interp_w)
            % Evaluates the oned dimensional integral  
            %   f = INTEGRAL(poly, interp_w)
            %
            %   interp_w    - Either the nodal values (piecewise) or 
            %                 corefficents (spectral)
            
            f = obj.int_W*interp_w;
        end
        
    end
    
end