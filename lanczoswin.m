function omega = lanczoswin(window_lenght, L)
%LANCZOSWIN Creates a Lanczos window
omega = zeros(1, window_lenght);
for i = 1 : window_lenght
    % Am memorat bucata care se repeta din formula pentru optimizare. Plus
    % avem un caz particular, cand avem o limita speciala de calculat.
    calc_aux = 2 * pi * (2 * (i - 1) - window_lenght + 1) / 2 / (window_lenght - 1);
    if calc_aux == 0
        omega(i) = 1; 
        % Aici mi-a dat batai de cap pana sa ma prin ca aveam
        % lim(x->0)sin(x)/x ...
    else
        % Problema de instabilitate numerica, pentru L subunitar, nu am
        % voie sa am numere complexe si trebuie tratat cazul. Nu ar trebui
        % sa aiba influenta asupra rezultatului final...
%         if calc_aux < 0 && L < 1
%             calc_aux = abs(calc_aux);
%         end
        % Se pare ca nu a mers cu varianta de mai devreme, asa ca o sa fac
        % elementul sa fie partea lui reala. Nu ar trebui sa fie gresit, mai ales
        % pentru acele valori mici. Din cate am observat pe datele
        % generate, partea imaginara e extrem de mica.
        omega(i) = real((sin(calc_aux) / calc_aux) ^ L);
    end
end
end