classdef pre_sol < full_sol
    properties
        precond
    end


    methods
        function sol = pre_sol(model, poly, opt, lowopt, N, epd, precond)
            sol@full_sol(model, 1, poly, opt, lowopt, N, epd);
            sol.precond = precond;
            sol.precond.SIRTs = cell(model.T, 1);
            sol.precond.L = sol.L;
            sol.precond.mu = sol.mu;
        end

        function sol = solve(sol)
            sol.samples(1:sol.model.d+sol.model.m,:,1) = priorsam(sol.model, sol.N);
            for t = 1:sol.model.T
                [temp, w_full] = push_samples(sol.model, sol.samples(1:sol.model.d+sol.model.m, :, t), ones(1, sol.N)/sol.N, t);
                sol.samples(:,:,t+1) = [temp;...
                    sol.samples(sol.model.d+1:sol.model.d+sol.model.m, :, t)];
                
                fprintf("FTT approximation at time %d because ESS =  %.2f / %d \n", t, ESS(w_full), sol.N)

                [sol, w] = reapprox(sol, w_full, t);

                % compute ESS
                sol.ESS_all(t) = ESS(w);
                sol.weight(:, t) = w';
            end
        end

        function [sol, w_full] = reapprox(sol, w_full, t)
            tic
            % find preconditioning parameters
            stepess = ESS(w_full);
            sam_new = sol.samples(:,:,t+1);
            w2 = w_full;
            ndouble = 0;

            if strcmp(sol.precond.opt, "g") || strcmp(sol.precond.opt, "geta")
                while stepess < sol.N/10 || ndouble == 0
                    if t == 1
                        ndouble = ndouble+1;
                        sam_old2 = priorsam(sol.model, 2^ndouble * sol.N);
                        [temp, ~] = push_samples(sol.model, sam_old2(1:sol.model.d+sol.model.m, :), ones(1, 2^ndouble * sol.N), t);
                        sam_new = [temp; sam_old2(sol.model.d+1:sol.model.d+sol.model.m, :)];
                        w_half = sol.precond.c*log(like(sol.model, sam_new(1:sol.model.d+sol.model.m,:), t));
                        w_half = w_half - max(w_half);
                        w_half = exp(w_half);
                        w_half(isnan(w_half)) = 0;
                        w2 = w_half/sum(w_half);
                        stepess = sum(w2)^2/sum(w2.^2);
                        if ndouble > 4
                            break
                        end
                    else
                        Tu2x = @(u) sol.precond.L(:, :, t-1)*eval_irt(sol.precond.SIRTs{t-1}, sol.precond.R(u))+sol.precond.mu(:, t-1);
                        ndouble = ndouble+1;
                        sam_old2 = sol.L(:,:,t-1) * eval_irt(sol.SIRTs{t-1}, rand(sol.model.d+2*sol.model.m, 2^ndouble * sol.N)) + sol.mu(:, t-1);
                        sam_old2 = Tu2x(sam_old2);
                        [temp, ~] = push_samples(sol.model, sam_old2(1:sol.model.d+sol.model.m, :), ones(1, 2^ndouble * sol.N), t);
                        sam_new = [temp; sam_old2(sol.model.d+1:sol.model.d+sol.model.m, :)];
                        w_half = sol.precond.c*log(like(sol.model, sam_new(1:sol.model.d+sol.model.m,:), t));
                        w_half = w_half - max(w_half);
                        w_half = exp(w_half);
                        w_half(isnan(w_half)) = 0;
                        w2 = w_half/sum(w_half);
                        stepess = sum(w2)^2/sum(w2.^2);
                        if ndouble > 4
                            break
                        end
                    end
                end
            else
                while stepess < sol.N/10
                    if t == 1
                        ndouble = ndouble+1;
                        sam_old2 = priorsam(sol.model, 2^ndouble * sol.N);
                        [temp, w2] = push_samples(sol.model, sam_old2(1:sol.model.d+sol.model.m, :), ones(1, 2^ndouble * sol.N), t);
                        sam_new = [temp; sam_old2(sol.model.d+1:sol.model.d+sol.model.m, :)];
                        stepess = sum(w2)^2/sum(w2.^2);
                        if ndouble > 4
                            break
                        end
                    else
                        Tu2x = @(u) sol.precond.L(:, :, t-1)*eval_irt(sol.precond.SIRTs{t-1}, sol.precond.R(u))+sol.precond.mu(:, t-1);
                        ndouble = ndouble+1;
                        sam_old2 = sol.L(:,:,t-1) * eval_irt(sol.SIRTs{t-1}, rand(sol.model.d+2*sol.model.m, 2^ndouble * sol.N)) + sol.mu(:, t-1);
                        sam_old2 = Tu2x(sam_old2);
                        [temp, w2] = push_samples(sol.model, sam_old2(1:sol.model.d+sol.model.m, :), ones(1, 2^ndouble * sol.N), t);
                        sam_new = [temp; sam_old2(sol.model.d+1:sol.model.d+sol.model.m, :)];
                        stepess = sum(w2)^2/sum(w2.^2);
                        if ndouble > 4
                            break
                        end
                    end
                end
            end
            fprintf('Enhanced ESS is %.2f. \n',  stepess)
            
            [mup_temp, Lp_temp, ~] = computeL(sam_new, w2);
            sam_new = datasample(sam_new', size(sam_new, 2), 'Weights', w2);
            sam_new = sam_new';
            Lp_temp = sol.epd * Lp_temp;
            sol.precond.L(:, :, t) =  Lp_temp;
            sol.precond.mu(:, t) = mup_temp;
            Lpre = @(x) Lp_temp*x + mup_temp;
            
            % give prior
            if t == 1
                logfun_prior = @(x) -log(priorpdf(sol.model, x));
            else
                sirt_temp = sol.SIRTs{t-1};
                sirt_pre = sol.precond.SIRTs{t-1};
                if sirt_temp.int_dir ~= 1, sirt_temp = marginalise(sirt_temp, 1); end
                if sirt_pre.int_dir ~= 1, sirt_pre = marginalise(sirt_pre, 1); end

                Lpreinv = @(x) sol.precond.L(1:sol.model.d+sol.model.m, 1:sol.model.d+sol.model.m, t-1)\...
                    (x-sol.precond.mu(1:sol.model.d + sol.model.m, t-1));
                Linv = @(x) sol.L(1:sol.model.d+sol.model.m, 1:sol.model.d+sol.model.m, t-1)\...
                    (x-sol.mu(1:sol.model.d + sol.model.m, t-1));
                Tx2u = @(x) sol.precond.Rinv(eval_rt(sirt_pre, Lpreinv(x))); 
                pre_pdfapp = @(x) eval_pdf(sirt_pre, Lpreinv(x));           
                fun_prior = @(x) eval_pdf(sirt_temp, Linv(Tx2u(x))).*pre_pdfapp(x)./ sol.precond.r(Tx2u(x));
                logfun_prior = @(x) -log(fun_prior(x));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% select the preconditioner and the corresponding preconditioned posterior

            switch sol.precond.opt
                case "pifg"
                    % These are (pi/eta)^c * eta
                    fun_pre = @(x, const) fun_post(sol, logfun_prior, t, Lpre(x), const) * sol.precond.c - ...
                        (1-sol.precond.c) * sol.precond.logr(x);
                    fun_pred = @(u, const, Tu2x, L_temp, mu_temp) ...
                        fun_post(sol, logfun_prior, t, Tu2x(L_temp*u + mu_temp), const) * (1-sol.precond.c) ...
                        +(1-sol.precond.c)* sol.precond.logr(Lp_temp\(Tu2x(L_temp*u + mu_temp)-mup_temp)) ...
                        - sol.precond.logr(L_temp*u + mu_temp);

                case "fg"
                    % These are pi_{k-1}(fg/eta)^c * eta
                    fun_pre = @(x, const) logfun_prior(x([1:sol.model.d, sol.model.d + sol.model.m+1:end], :))...
                        - log(transition(sol.model, x, t))* sol.precond.c - log(like(sol.model, x, t)) * sol.precond.c- const;

                    fun_pre = @(x, const) fun_pre(Lpre(x), const) - (1-sol.precond.c) * sol.precond.logr(x) ;
                    fun_pred = @(u, const, Tu2x, L_temp, mu_temp) ...
                        - log(transition(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c)...
                        - log(like(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c) -const ...
                        +(1-sol.precond.c)* sol.precond.logr(Lp_temp\(Tu2x(L_temp*u + mu_temp)-mup_temp)) ...
                        - sol.precond.logr(L_temp*u + mu_temp);

                case "fgeta"
                    fun_pre = @(x, const) logfun_prior(x([1:sol.model.d, sol.model.d + sol.model.m+1:end], :))...
                        - log(transition(sol.model, x, t))* sol.precond.c - log(like(sol.model, x, t)) * sol.precond.c- const;

                    fun_pre = @(x, const) fun_pre(Lpre(x), const);
                    fun_pred = @(u, const, Tu2x, L_temp, mu_temp) ...
                        - log(transition(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c)...
                        - log(like(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c) -const ...
                        - sol.precond.logr(L_temp*u + mu_temp);

                case "g"
                    % These are pi_{k-1}f(g/eta)^c * eta
                    fun_pre = @(x, const) logfun_prior(x([1:sol.model.d, sol.model.d + sol.model.m+1:end], :))...
                        - log(transition(sol.model, x, t)) - log(like(sol.model, x, t)) * sol.precond.c- const;

                    fun_pre = @(x, const) fun_pre(Lpre(x), const);
                    fun_pred = @(u, const, Tu2x, L_temp, mu_temp) ...
                        - log(like(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c) -const ...
                        - sol.precond.logr(L_temp*u + mu_temp);
                    
                case "geta"
                    fun_pre = @(x, const) logfun_prior(x([1:sol.model.d, sol.model.d + sol.model.m+1:end], :))...
                        - log(transition(sol.model, x, t)) - log(like(sol.model, x, t)) * sol.precond.c- const;

                    fun_pre = @(x, const) fun_pre(Lpre(x), const) - (1-sol.precond.c) * sol.precond.logr(x) ;
                    fun_pred = @(u, const, Tu2x, L_temp, mu_temp) ...
                        - log(like(sol.model, Tu2x(L_temp*u + mu_temp), t)) * (1-sol.precond.c) -const ...
                        +(1-sol.precond.c)* sol.precond.logr(Lp_temp\(Tu2x(L_temp*u + mu_temp)-mup_temp)) ...
                        - sol.precond.logr(L_temp*u + mu_temp);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % establish the preconditioner
            bound = sol.poly.oned_domains{1}.bound(2);
            if isa(sol.poly.oned_domains{1}, "AlgebraicMapping")
                bound = sol.poly.oned_domains{1}.scale;
            end
            z = rand(sol.model.d + 2*sol.model.m, 4e3) *2 * bound - bound;
            fun_const = @(x) nan2inf(fun_pre(x, 0));
            const = min(fun_const(z));
            logfun_pre = @(x) nan2inf(fun_pre(x, const));

            % resample unweighted samples to serve as sample points and
            % debugger in TTSIRT
            samples_unweighted = Lp_temp\(sam_new - mup_temp);
            samples_debug = samples_unweighted(:, 1:size(sam_new,2)/2);
            samples_init = samples_unweighted(:, size(sam_new,2)/2+1:end);
            deb = InputData(samples_init, samples_debug);

            if t == 1
                sirt = TTSIRT(logfun_pre, sol.poly, sol.opt, 'var', deb);
            else
                sirt_temp = sol.precond.SIRTs{t-1};
                sirt = TTSIRT(logfun_pre,  sirt_temp, 'var', deb);
            end

            sol.precond.SIRTs{t} = sirt;
            Tu2x = @(u) Lpre(eval_irt(sirt, sol.precond.R(u)));
            Tx2u = @(x) sol.precond.Rinv(eval_rt(sirt, Lp_temp\(x-mup_temp))); 

            L_temp = eye(sol.model.d+2*sol.model.m);
            mu_temp = zeros(sol.model.d+2*sol.model.m, 1);
            
            sol.L(:, :, t) =  L_temp;
            sol.mu(:, t) = mu_temp;


            % establish fun_into_sirt
            bound = sol.poly.oned_domains{1}.bound(2);
            if isa(sol.poly.oned_domains{1}, "AlgebraicMapping")
                bound = sol.poly.oned_domains{1}.scale;
            end
            z = rand(sol.model.d + 2*sol.model.m, 4e3) *2 * bound - bound;
            fun_const = @(u) nan2inf(fun_pred(u, 0, Tu2x, L_temp, mu_temp));
            
            zz = fun_const(z);
            const = min(zz);
            logfun_post = @(u) nan2inf(fun_pred(u, const, Tu2x, L_temp, mu_temp));

            % resample unweighted samples to serve as sample points and
            % debugger in TTSIRT
            samples_unweighted = L_temp\(Tx2u(sam_new)- mu_temp);                
            samples_debug = samples_unweighted(:, 1:size(sam_new,2)/2);
            samples_init = samples_unweighted(:, size(sam_new,2)/2+1:end);
            deb = InputData(samples_init, samples_debug);

            sirt_temp = sol.precond.SIRTs{t};
            sirt_temp.approx.opt = sol.lowopt;
            sirt = TTSIRT(logfun_post,  sirt_temp, 'var', deb);

            % generate samples
            z = rand(sol.model.d+2*sol.model.m, sol.N);
            [r, p] = eval_irt(sirt, z);
            r_temp = sol.L(:,:,t) * r + sol.mu(:, t);
            sol.samples(:, :, t+1) = Tu2x(r_temp);

            pdf_true = fun_post(sol, logfun_prior, t, sol.samples(:, :, t+1), 0);
            
            pdf_app = p - log(eval_pdf(sol.precond.SIRTs{t}, ...
                Lp_temp\(sol.samples(:, :, t+1)-mup_temp))) + sol.precond.logr(r_temp);
            w_full = exp(pdf_app - pdf_true);
            fprintf("The total ESS is %.2f out of %d. \n", ESS(w_full), length(w_full))

            sol.FTT_time = [sol.FTT_time ,toc]; 
            % disp(sol.FTT_time)
            sol.SIRTs{t} = sirt;
        end

        function logpdf = fun_post(sol, logf, t, x, const)
            logpdf = logf(x([1:sol.model.d, sol.model.d + sol.model.m+1:end], :))...
                - log(transition(sol.model, x, t)) - log(like(sol.model, x, t)) - const;
                logpdf(isnan(logpdf)) = inf;
            % end
        end

        function [thetas, sams, w, varargout] = smooth(sol, N, T)
            d = sol.model.d;
            m = sol.model.m;
            sams = zeros(m , N, T+1);
            logpdf_eall = zeros(N, T);

            sirt = sol.SIRTs{T};
            [r,f] = eval_irt(sirt, rand(d+2*m, N));
            f = f - min(f);
            r = sol.L(:,:,T) *r + sol.mu(:, T);
            temp = eval_irt(sol.precond.SIRTs{T}, sol.precond.R(r));
            temp2 = sol.precond.L(:, :, T) * temp + sol.precond.mu(:, T); % recover samples using L and mu
            nlogpdf_e = f - log(eval_pdf(sol.precond.SIRTs{T}, temp)) + sol.precond.logr(r);

            sams(:, :, T+1) = temp2(d+1: d+m, :); % generate smaples
            sams(:, :, T) = temp2(d+m+1: end, :); % generate smaples
            if d > 0
                thetas = temp2(1:d, :);
            else
                thetas = [];
            end

            logpdf_e = -nlogpdf_e;
            logpdf_eall(:, T) = logpdf_e;
            for t = T-2:(-1):0  % k is time, stored in k+1(second indice) in samples
                disp(t)
                sirt = sol.SIRTs{t+1};
                sirtp = sol.precond.SIRTs{t+1};
                xt = sams(:, :, t+2);
                L = @(x) sol.L(:,:, t+1)*x + sol.mu(:,t+1);
                Linv = @(x) sol.L(1:d+m,1:d+m, t+1)\(x - sol.mu(1:d+m,t+1));
                Lpre = @(x) sol.precond.L(:,:, t+1)*x + sol.precond.mu(:,t+1);
                Lpreinv = @(x) sol.precond.L(1:d+m,1:d+m, t+1)\(x - sol.precond.mu(1:d+m,t+1));
                if d > 0
                    z = Linv(sol.precond.Rinv(eval_rt(sirtp, Lpreinv([thetas; xt]))));
                else
                    z = Linv(sol.precond.Rinv(eval_rt(sirtp, Lpreinv(xt))));
                end
                [r,f] = eval_cirt(sirt, z, rand(m, N));
                f = f - min(f);

                temp = L([z ; r]);
                temp2 = eval_irt(sirtp, sol.precond.R(temp));
                temp3 = Lpre(temp2);
                sams(:, :, t+1) = temp3(d+m+1:end, :);
                logpdf_e = logpdf_e - f + log(eval_pdf(sirtp, temp2)) - log(eval_pdf(sirtp, temp2(1:d+m, :)))...
                    + sol.precond.logr(temp(1:d+m, :)) - sol.precond.logr(temp);
                logpdf_eall(:, t+1) = logpdf_e;
            end

            
            % theoretical pdf

            fulldata = zeros(d+2*m, N, T+1);
            if d > 0
                fulldata(1:d, :, :) = repmat(thetas, 1, 1, T+1);
            end
            fulldata(d+1:d+m, :, :) = sams;
            fulldata(d+m+1:end, :, 2:T+1) = sams(:,:, 1:T);

            logpdf_t = log(priorpdf(sol.model, fulldata(1:d+m, :, 1)));
            for k = 1:T
                logpdf_t = logpdf_t + log(transition(sol.model, fulldata(:, :, k+1), k)) +...
                    log(like(sol.model, fulldata(:, :, k+1), k));
            end
            
            % compute ESS
            w = logpdf_t - logpdf_e;
            w = exp(w - max(w));
            w = w/sum(w);
            if nargout == 4
                varargout{1} = logpdf_eall;
            end
        end

        function ESSs = smooth_t(sol, thetas, sams, logpdf_eall, Tts)
            d = sol.model.d;
            m = sol.model.m;
            T = size(sams, 3) - 1;
            N = size(sams, 2);
            ESSs = zeros(size(Tts));

            for k = 1:length(Tts)
                t_init = Tts(k);
                logpdf_e = logpdf_eall(:, t_init)';
                % theoretical pdf

                fulldata = zeros(d+2*m, N, T-t_init+1);
                if d > 0
                    fulldata(1:d, :, :) = repmat(thetas, 1, 1, T-t_init+1);
                end
                fulldata(d+1:d+m, :, :) = sams(:, :, t_init+1:T+1);
                fulldata(d+m+1:end, :, :) = sams(:,:, t_init:T);

                temp = sol.precond.L(:,:, t_init)\(fulldata(:, :, 1)-sol.precond.mu(:, t_init));
                r = sol.precond.Rinv(eval_rt(sol.precond.SIRTs{t_init}, temp));
                r2 = sol.L(:,:, t_init)\(r-sol.mu(:, t_init));
                f = log(eval_pdf(sol.SIRTs{t_init}, r2));
                logpdf_t = f + log(eval_pdf(sol.precond.SIRTs{t_init}, temp)) - sol.precond.logr(r);

                for t = 1:T-t_init
                    logpdf_t = logpdf_t + log(transition(sol.model, fulldata(:, :, t+1), t_init+t)) +...
                        log(like(sol.model, fulldata(:, :, t+1), t_init+t));
                end

                w = logpdf_t - logpdf_e;
                w = exp(w - max(w));
                w = w/sum(w);
                ESSs(k) = ESS(w);
            end
        end

        function plot_sirt(sol, t)
            d = sol.model.d;
            m = sol.model.m;
            % sams = zeros(d + 2*m , N);

            sirt = sol.SIRTs{t};
            sirtp = sol.precond.SIRTs{t};
            [r,~] = eval_irt(sirt, rand(d+2*m, 1e3));
            sams = sol.L(:, :, t) * r + sol.mu(:, t); % recover samples using L and mu
            sams = eval_irt(sirtp, sol.precond.R(sams));
            sams = sol.precond.L(:, :, t) * sams + sol.precond.mu(:, t);
            
            plot_marginals(compute_marginals(sams))
        end

        function stats = plot_stats(sol, thetas, sams, w)
            stats = plot_stats@full_sol(sol, thetas, sams, w);
        end
    end
end