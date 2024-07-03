classdef red_sol < Y_sol
    properties
        sqr % indicate if squared version of FTT is used
        cov_sigma % store the transformed cov of \bEp
        lowopt
    end


    methods
        function sol = red_sol(model, sqr, poly, opt, lowopt, N, epd)
            sol@Y_sol(model, poly, opt, N, epd, 10);
            sol.lowopt = lowopt;
            sol.sqr = sqr;
            sol.SIRTs = cell(model.T, 1);
            sol.cov_sigma = sol.model.pre.errcov;
        end

        function sol = solve(sol)
            sol = solve@Y_sol(sol);
        end

        function [sol, logfun_post] = reapprox(sol, w, t, ~)
            % convolution
            tic
            samples_unweighted = datasample(sol.samples(:,:,t+1),...
                sol.N, 2, 'Weights', w);

            mu = sum(sol.samples(:,:,t+1).*w, 2);
            cov_x = sol.N/(sol.N-1) * (w.*(sol.samples(:,:,t+1) - mu))*(sol.samples(:,:,t+1) - mu)';
            
            
            if ESS(w) < 0.1*sol.N
                [L, cov_k] = precond(1.5*sol.epd^2*cov_x, sol.model.pre.errcov);
            else
                [L, cov_k] = precond(sol.epd^2*cov_x, sol.model.pre.errcov);
            end
            sol.L(:, :, t) = L;
            sol.mu(:, t) = mu;
            
            
            % give prior
            if isempty(which('update_cov'))
                fun_updatecov = @(x) diag(sol.cov_sigma).*ones(1, size(x, 2));
            else
                fun_updatecov = @(x) update_cov(sol.cov_sigma, x, sol.L(:,:,t-1), sol.mu(:, t-1));
            end

            if t == 1
                logfun_prior = @(x) -log(priorpdf_conv(sol.model, x));
            else
                if sol.sqr == 1
                    [cores_new, poly_new] = TT_expand(sol.SIRTs{t-1});
                    logfun_prior = @(x) -log(TT_eval(cores_new, poly_new,...
                        fun_updatecov, sol.L(:, :, t-1)\(x- sol.mu(:, t-1))));
                else
                    ftt = sol.SIRTs{t-1}.approx; % this is a TTFun, with fields base and data
                    logfun_prior = @(x) -log(max(TT_conv_Lag1(ftt, fun_updatecov, sol.L(:, :, t-1)\(x- sol.mu(:, t-1))), 0));
                end
            end

            % establish fun_into_sirt
            const = 0;
            bound = sol.poly.oned_domains{1}.bound(2);
            z = rand(sol.model.d + sol.model.m, 1e3) *2 * bound - bound;
            fun_const = @(x) fun_into_sirt(sol, logfun_prior, t, t-1, L*x + mu, const);
            const = min(fun_const(z));
            logfun_post = @(x) fun_into_sirt(sol, logfun_prior, t, t-1, L*x + mu, const);

            samples_unweighted = L\(samples_unweighted - mu);
            samples_debug = samples_unweighted(:, 1:sol.N/2);
            samples_init = samples_unweighted(:, sol.N/2+1:end);
            deb = InputData(samples_init, samples_debug);

            % build IRT or SIRT
            
            if sol.sqr == 0
                if t == 1
                    ftt = TTIRT(logfun_post, sol.poly, sol.opt, 'var', deb);
                else
                    ftt_temp = sol.SIRTs{t-1};
                    ftt_temp.approx.opt = sol.lowopt;
                    ftt = TTIRT(logfun_post, ftt_temp, 'var', deb);
                end
                sol.SIRTs{t} = ftt;
            else
                if t == 1
                    sirt = TTSIRT(logfun_post, sol.poly, sol.opt, 'var', deb);
                else
                    sirt_temp = sol.SIRTs{t-1};
                    sirt_temp.approx.opt = sol.lowopt;
                    sirt = TTSIRT(logfun_post,  sirt_temp, 'var', deb);
                end
                sol.SIRTs{t} = sirt;
            end

            sol.cov_sigma = cov_k;
            sol.FTT_time = [sol.FTT_time ,toc]; 
        end

        function logpdf = fun_into_sirt(sol, logf, t, T_r, x, const)
            logpdf = fun_into_sirt@Y_sol(sol, logf, t, T_r, x, const);
        end

    end
end