clc; close all; clear;

tic

%% Definice velikosti rastru
rx = 50;
ry = rx;

%% Definice testovacích krivek
T = 2*pi;
% t = 0:5e-3:T;
t = linspace(0, T, 3000);
lx = cos(t); ly = sin(t); lz = abs(sin(3*t/2));
% lx = 3*cos(t) + cos(3*t); ly = 3*sin(t) + sin(3*t); lz = ones(1,size(t, 2)); % nephroid
% lx = cos(t); ly = sin(t).*cos(t); lz = ones(1,size(t, 2)); % infinity
% jj=3; kk=5; lx = cos(t) - cos(t).^jj; ly = sin(t) - sin(t).^kk; lz = ones(1,size(t, 2)); %3 3, 3 5

[b, maskbound, mask, scalepar] = getBound(lx, ly, lz, rx, ry); %% zobrazení pocátecních podmínek do pracovního prostoru
A = getA(ry, rx); % Generování matice A
S = reshape(b, [], 1); % Serializace pracovního prostoru
 
maxe = 1e-11; % nejmenší progres, u kterého se algoritmus zastaví
err = 1; % první progress, jen pro podmínku ve while
% p = [160, 100];
p = round([rx/2, ry/2]); % bod v rastru, kde se merí progress
pnum = round(p(1)*rx + p(2)); % serializace sou?adnic bodu
tmp = 5; % pro první kontrolu while podmínky
iter = 0;

% figure
% pause(5);
while (abs(S(pnum) - tmp) > maxe || S(pnum) < 0.1) %% zastavení algoritmu pri progresu menším než "maxe"
    % while(iter < 3000)
    tmp = S(pnum);
    S = A*S;
    S(maskbound) = b(maskbound);
    
    iter = iter + 1;
%     if (mod(iter, 100) == 0) %% pro postupné vykreslování a ukládání progresu
%         err = [err, (S(pnum) - tmp)];
%         K = rot90(rot90(reshape(S, ry, rx)));
%         K(~mask) = NaN;
%         hh = surf(K, 'EdgeColor', 'none');
%         title(['Minimal surface - Jacobi method, ' int2str(iter) ...
%             ' iterations, with time: '...
%             num2str(toc) 's, delta progress in [0,0]= ' num2str((S(pnum) - tmp))], 'FontSize', 15);
%         hh.XData = linspace(min(lx), max(lx), rx);
%         hh.YData = linspace(min(ly), max(ly), ry);
%         xlabel('x'); ylabel('y'); zlabel('z')
%         drawnow;
%     end
end

endtime = toc

S = reshape(S, ry, rx);
S(~mask) = NaN;
hh = surf(rot90(rot90(S)), 'EdgeColor', 'none'); %view(0, 90);
hh.XData = linspace(min(lx), max(lx), rx);
hh.YData = linspace(min(ly), max(ly), ry);

axmin = min(min(lx), min(ly));
axmax = max(max(lx), max(ly));
axis([axmin, axmax, axmin, axmax, min(min(S)), max(max(S))]);

hold on
plot3(max(lx)*lx/(scalepar), max(ly)*ly/(scalepar), lz, 'LineWidth', 5, 'Color', 'red');
hold off
title(['Jacobi method, time: ' num2str(endtime*1e3) ' ms, iterations: ' num2str(iter)], 'FontSize', 16);
xlabel('x', 'FontSize', 14); ylabel('y', 'FontSize', 14); zlabel('z', 'FontSize', 14)

% figure; plot(linspace(0, iter, numel(err(2:end))), err(2:end), 'LineWidth', 3); 
% title('Delta progress function', 'FontSize', 16);
% xlabel('Iterations [-]', 'FontSize', 14); ylabel('Diference [-]', 'FontSize', 14);
% grid on

