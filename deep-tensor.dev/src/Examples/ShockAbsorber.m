classdef ShockAbsorber

    properties
        y
        censind
        D
        x
        beta_mean
        beta_var
        ref
    end

    methods
        function obj = ShockAbsorber(varargin)

            param = struct;
            for i=1:2:numel(varargin)
                param.(varargin{i}) = varargin{i+1};
            end

            % Data table from the paper
            obj.y = [6700 6950 7820 8790 9120 9660 9820 11310 11690 11850 11880 12140 12200 ...
                12870 13150 13330 13470 14040 14300 17520 17540 17890 18420 18960  ...
                18980 19410 20100 20100 20150 20320 20900 22700 23490 26510 27410  ...
                27490 27890 28100];
            obj.censind = [0 1 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 0 1 0 1 0 1 1];

            if (~isfield(param, 'D'))
                obj.D = input('Number of covariates D = ? (default 6): ');
                if (isempty(obj.D))
                    obj.D = 6;
                end
            else
                obj.D = param.D;
            end
            if (~isfield(param, 'x'))
                obj.x = input('Covariates data x = ? (default empty, means random generation): ');
            else
                obj.x = param.x;
            end

            if (isempty(obj.x))
                % Simulate artificial covariates
                obj.x = randn(obj.D, numel(obj.y))*1/obj.D;
                fprintf('\tcovariates: randomly synthesized, mean(x)=%g\n', mean(obj.x(:)));
            else
                fprintf('\tcovariates: given, mean(x)=%g\n', mean(obj.x(:)));
            end

            % scales for beta_0
            obj.beta_mean = log(30796);
            obj.beta_var = 0.1563;

            % Gaussian reference
            obj.ref = GaussReference(0, 1, UnboundedDomain(0,1));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkd, mlp] = minus_log_post(obj, theta)
            theta = theta';
            theta(:,end) = theta(:,end)+1;
            mllkd = shock_log_weibull(obj, theta);
            mllkd = -mllkd';
            mlp = shock_log_prior(obj, theta);
            mlp = -mlp';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Weibull likelihood
        function F = shock_log_weibull(obj, theta)
            d = size(theta,2)-2;
            beta = theta(:,1:d+1);
            lambda = exp(theta(:,d+2));
            F = zeros(size(theta,1), 1);
            m = numel(obj.y);
            for i=1:m
                logeta = beta(:,1)*obj.beta_var + obj.beta_mean + beta(:,2:d+1)*obj.x(:,i);
                eta = exp(logeta);
                yeta = obj.y(i)./eta;
                if (obj.censind(i))
                    f = -(yeta.^lambda); % censored PDF = 1-CDF
                else
                    f = log(lambda) - logeta + (lambda-1).*(log(obj.y(i)) - logeta) - (yeta.^lambda);
                    f = f + log(3e4); % to prevent underflow
                end
                F = F + f;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % prior distribution for shock absorber
        function F = shock_log_prior(obj, theta)
            d = size(theta,2)-2;
            % Normal-Gamma
            %alpha = 6.8757; % original setup
            %beta = 2.2932; % original setup
            %
            alpha = 0.5; % original setup
            beta = 1E-4; % original setup
            theta2 = log(exp(theta(:,end))+1);
            logJ = theta(:,end) - theta2;
            %F = (alpha-0.5+1)*theta2 - exp(theta2).*(beta+sum((theta(:,1:d+1).^2),2)/2);
            F = logJ + (alpha-0.5)*log(theta2) - theta2.*(beta+sum((theta(:,1:d+1).^2),2)/2);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function f = potential(obj, theta)
            [mllkd, mlp] = minus_log_post(obj, theta);
            f =  mllkd+mlp;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkd,mlp] = pullback_potential_pcn(obj, irt, z)
            u = eval_cdf(obj.ref, z);
            [x,mlogfx] = eval_irt(irt, u);
            [mllkdx, mlpx] = minus_log_post(obj, x);

            mllkd = mllkdx+mlpx-mlogfx;
            mlp = 0.5*sum(z.^2,1);
        end

    end
end