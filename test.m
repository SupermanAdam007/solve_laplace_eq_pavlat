clc; close all; clear;

statusm1 = [];
statusm2 = [];
% rastrpos = 1;

% figure; grid on
for rx = 10:10:600
    %     maxepos = 1;
    %     for maxe = 1e-4:-5e-6:1e-7
    ry = rx;
    T = 2*pi;
    t = linspace(0, T, 3000);
    lx = cos(t); ly = sin(t); lz = abs(sin(3*t/2));
    
    [b, maskbound, mask, scalepar] = getBound(lx, ly, lz, rx, ry); %% zobrazení pocátecních podmínek do pracovního prostoru
    
    tic
    A = getA(ry, rx); % Generování matice A
    S = reshape(b, [], 1); % Serializace pracovního prostoru
    
    maxe = 1e-2; % nejmenší progres, u kterého se algoritmus zastaví
    p = round([rx/2, ry/2]); % bod v rastru, kde se merí progress
    pnum = round(p(1)*rx + p(2)); % serializace sou?adnic bodu
    tmp = 5; % pro první kontrolu while podmínky
    while (abs(S(pnum) - tmp) > maxe || S(pnum) < 0.1) %% zastavení algoritmu pri progresu menším než "maxe"
        tmp = S(pnum);
        S = A*S;
        S(maskbound) = b(maskbound);
    end
    statusm1 = [statusm1 toc];
    %         status{rastrpos, maxepos} = {rastrpos, maxepos};
    %         statusm(rastrpos, maxepos) = toc;
    %         maxepos = maxepos + 1;
    
    
    %% konvoluce
    tic
    tmp = 5;
    lap = 1/4*[0 1 0; 1 0 1; 0 1 0];
    iter = 0; res = b;
    while ((abs(res(p(1), p(2)) - tmp) > maxe) | (res(p(1), p(2)) < 0.1))
        tmp = res(p(1), p(2));
        res = conv2(res, lap, 'same');
        res(maskbound) = b(maskbound);
    end
    statusm2 = [statusm2 toc];
    
    fprintf('Result1: rx = %d, time1 = %4.2f s, time2 = %4.2f s\n', rx, statusm1(end), statusm2(end));
%     plot(statusm1, statusm2); legend('Serialization', 'convolution'); drawnow;
    %     end
    %     rastrpos = rastrpos + 1;
end
