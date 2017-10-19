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
    
    [b, maskbound, mask, scalepar] = getBound(lx, ly, lz, rx, ry); %% zobrazen� poc�tecn�ch podm�nek do pracovn�ho prostoru
    
    tic
    A = getA(ry, rx); % Generov�n� matice A
    S = reshape(b, [], 1); % Serializace pracovn�ho prostoru
    
    maxe = 1e-2; % nejmen�� progres, u kter�ho se algoritmus zastav�
    p = round([rx/2, ry/2]); % bod v rastru, kde se mer� progress
    pnum = round(p(1)*rx + p(2)); % serializace sou?adnic bodu
    tmp = 5; % pro prvn� kontrolu while podm�nky
    while (abs(S(pnum) - tmp) > maxe || S(pnum) < 0.1) %% zastaven� algoritmu pri progresu men��m ne� "maxe"
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
