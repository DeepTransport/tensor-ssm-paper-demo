function dh2 = test_irt(func, irt, x)

if isscalar(x)
    x = random(irt, ceil(x));
end

mlogfx = eval_potential(irt,x);
mlpt = func(x);

[~,dh2] = f_divergence(-mlpt,-mlogfx);

end