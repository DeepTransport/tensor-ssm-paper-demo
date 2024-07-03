function pdf = like(model, thetax, t)

y = model.Y(t);
beta =  normcdf(thetax(2,:))*0.8+0.1;
cc = beta .* exp(0.5*thetax(3,:));
pdf = normpdf(y./cc)./cc;
end