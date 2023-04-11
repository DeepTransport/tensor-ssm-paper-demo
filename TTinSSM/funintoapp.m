function val = funintoapp(u, mu, L, That, rt2, pdfth, pdfthold, model, k)
% this is the function approximated in FTT in DIRT
% the input sirt is to compute pdf of x conditional on theta

txx = That(L*u+mu);
val = (transition(txx, model).*...
         like(txx, model, model.Y(:,k)).* ...
         pdfthold(rt2*txx) ./ model.fu(model.L_theta(:,:,k-1)\(txx - model.mu_theta(:,k-1)))).^(1-model.c).*model.fu(L*u+mu);
end