function [Delta_pr, Delta_sr] = check_PPFTI(h, omega_p, omega_s)
%CHECK_PPFTI Functie pentru rezolvarea PPFTI
%   Este calculat dupa modalitatea descrisa in fisierul de indicatii.
%   Parametrul res_freqz este doar ca sa setez rezolutia.
frecv_p = linspace(0, omega_p, 2000);
frecv_s = linspace(omega_s, pi, 2000);
H = freqz(h, 1, frecv_p);
Delta_pr = max(abs(1 - abs(H)));
H = freqz(h, 1, frecv_s);
Delta_sr = max(abs(H));
end

