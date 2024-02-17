function [h, omega_c, window_length, delta_pr, delta_sr] = best_filter_4(omega_p, omega_s, Delta_p, Delta_s, tip, parametru)
%BEST_FILTER_4 Genereaza cel mai bun filtru specific unei ferestre.
%   Folosind o rutina specifica pentru fiecare tip de fereastra pentru
%   generare (avand in vedere ca putem avea un parametru specific unui tip
%   de fereastra). Apoi acesta evalueaza care ar fi fereastra cea mai buna,
%   bazata pe lungimea ferestrei, apoi pe tolerantele pe care le are
%   fereastra.
%   Valorile pentru tip sunt tukey, lanczos, blackman, hamming, hanning,
%   chevysev, kaiser, triunghiular si dreptunghiular.
%   Este posibil pentru celelalte ferestre sa gasim rezultate
%   surprinzatoare (cand criteriul principal este dimensiunea ferestrei).
arguments
    omega_p (1, 1) double
    omega_s (1, 1) double
    Delta_p (1, 1) double
    Delta_s (1, 1) double
    tip string
    parametru (1, 1) double = 0
end
notfound = true;
window_length = 4;
if parametru < 0 
    error('Parametrul nu are cum sa fie negativ.')
else
    while notfound
        w_cut_weight = 0.5;
        omega_c = freq_cut_calc(omega_p, omega_s, w_cut_weight);
        freq_c = omega_c / pi;
        switch tip
            case 'dreptunghiular'
                w = rectwin(window_length);
            case 'triunghiular'
                w = triang(window_length);
            case 'blackman'
                w = blackman(window_length);
            case 'chebysev'
                w = chebwin(window_length, parametru);
%                 if parametru < 80 || parametru > 100
%                     warning('In afara valorilor uzuale')
%                 end
            case 'hamming'
                w = hamming(window_length);
            case 'hanning'
                w = hann(window_length);
            case 'kaiser'
                w = kaiser(window_length, parametru);
%                 if parametru < 2 || parametru > 10
%                     warning('In afara valorilor uzuale')
%                 end
            case 'lanczos'
                w = lanczoswin(window_length, parametru);
%                 if parametru < 0.8 || parametru > 1.2
%                     warning('In afara valorilor uzuale')
%                 end
            case 'tukey'
                if parametru > 100
                    error('Alpha nu are cum sa fie mai mare de 100%.')
                else
                    w = tukeywin(window_length, parametru);
                end
            otherwise
                error('Tipul de fereastra nu exista/nu a fost programat!')
        end
        h = fir1(window_length - 1, freq_c, w);
        [delta_pr, delta_sr] = check_PPFTI(h, omega_p, omega_s);
        if delta_pr > Delta_p || delta_sr > Delta_s
            window_length = window_length + 1;
        else
             % Urmatorul caz este pentru a incerca sa extragem fereastra cu frecventa
             % de taiere optima. Avand in vedere rezultatul generat la faza
             % 3, cel mai probabil este ca frecventa de taiere optima sa
             % fie tot in mijlocul zonei de tranzitie.
             best_weight = 0.5;
             sum_delta = delta_pr + delta_sr;
             for wei = 0.01 : 0.01 : 0.99
                 om_c = freq_cut_calc(omega_p, omega_s, wei);
                 freq_c = om_c / pi;
                 h = fir1(window_length - 1, freq_c, w);
                 [delta_pr, delta_sr] = check_PPFTI(h, omega_p, omega_s);
                 if delta_pr < Delta_p && delta_sr < Delta_s && sum_delta > (delta_pr + delta_sr)
                     best_weight = wei;
                     sum_delta = delta_sr + delta_pr;
                 end
             end
             notfound = false;
        end
    end
    omega_c = freq_cut_calc(omega_p, omega_s, best_weight); 
    h = fir1(window_length - 1, freq_c, w); 
    [delta_pr, delta_sr] = check_PPFTI(h, omega_p, omega_s);
end
end