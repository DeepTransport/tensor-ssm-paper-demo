function sol = rom_solve(rom, vr)

if rom.sq_param
    switch rom.kappa_type
        case{'scalar'}
            c = rom.Ks{1}*exp(rom.redu_us{1}*vr + rom.redu_means{1});
            tmp = c*c';
            x = tmp(rom.inds{1});
            Ak = reshape(rom.As{1}*x, rom.dof, rom.dof);
        case{'vector'}
            Ak = zeros(rom.dof);
            for di = 1:length(rom.As)
                ind_di = rom.redu_param_ind{di};
                c = rom.Ks{di}*exp(rom.redu_us{di}*vr(ind_di) + rom.redu_means{di});
                tmp = c*c';
                x = tmp(rom.inds{di});
                Ak = Ak + reshape(rom.As{di}*x, rom.dof, rom.dof);
            end
        case{'tensor'}
            c = rom.Ks{1}*exp(rom.redu_us{1}*vr(rom.redu_param_ind{1}) + rom.redu_means{1});
            tmp = c*c';
            x = tmp(rom.inds{1});
            Ak = reshape(rom.As{1}*x, rom.dof, rom.dof);
            % cross
            c = vr(rom.redu_param_ind{2});
            tmp = c*c';
            x = tmp(rom.inds{2});
            Ak = Ak + reshape(rom.As{2}*x, rom.dof, rom.dof);
    end
else
    switch rom.kappa_type
        case{'scalar'}
            x = rom.Ks{1}*exp(rom.redu_us{1}*vr + rom.redu_means{1});
            Ak = reshape(rom.As{1}*x, rom.dof, rom.dof);
        case{'vector'}
            Ak = zeros(rom.dof);
            for di = 1:length(rom.As)
                ind_di = rom.redu_param_ind{di};
                x = rom.Ks{di}*exp(rom.redu_us{di}*vr(ind_di) + rom.redu_means{di});
                Ak = Ak + reshape(rom.As{di}*x, rom.dof, rom.dof);
            end
    end
end

sol.state = (rom.A + Ak)\rom.b;
sol.d  = rom.obs_operator*sol.state;

% apply qoi function
if rom.qoi_flag
    sol.qoi = rom.phi*(Ak*sol.state);
else
    sol.qoi = [];
end

end