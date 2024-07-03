function ess = ESS(w)
    ess = sum(w, 'omitnan')^2/sum(w.^2, 'omitnan');
end