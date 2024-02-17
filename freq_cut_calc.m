function omega_c = freq_cut_calc(omega_p, omega_s, weight)
%FREQ_CUT_CALC Calculeaza frecventa de taiere
% Se ia dupa un parametru de pondere.
if weight < 0 || weight > 1
    error('Ponderea trebuie sa fie in [0, 1].')
else
    omega_c = (1 - weight) * omega_p + weight * omega_s;
end
end

