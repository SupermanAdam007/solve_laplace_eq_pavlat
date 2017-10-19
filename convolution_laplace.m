clc; close all; clear;

tic
%% Definice velikosti rastru
%         tic
rx = 50;
ry = rx;

%% Definice testovacích funkcí
T = 2*pi;
% t = 0:5e-3:T;
t = linspace(0, T, 1000);
lx = cos(t); ly = sin(t); lz = abs(sin(3*t/2));
% lx = 3*cos(t) + cos(3*t); ly = 3*sin(t) + sin(3*t); lz = ones(1,size(t, 2)); % nephroid
% lx = cos(t); ly = sin(t).*cos(t); lz = ones(1,size(t, 2)); % infinity
% jj=3; kk=5; lx = cos(t) - cos(t).^jj; ly = sin(t) - sin(t).^kk; lz = ones(1,size(t, 2)); %3 3, 3 5

[b, maskbound, mask, scalepar] = getBound(lx, ly, lz, rx, ry);

%%
res = b;
% figure; surf(res, 'EdgeColor', 'none'); view(-70, 40); axis([0 rx 0 ry 0 max(res(:))])

lap = 1/4*[0 1 0;
    1 0 1;
    0 1 0];

maxe = 1e-11; err = 1;
% p = [160 100];
% p = [140 150];
p = round([rx/2, ry/2]);
tmp = 5; iter = 0;
while ((abs(res(p(1), p(2)) - tmp) > maxe) | (res(p(1), p(2)) < 0.1))
    tmp = res(p(1), p(2));
    res = conv2(res, lap, 'same');
    res(maskbound) = b(maskbound);
    
%     iter = iter + 1;
%     if (mod(iter, 100) == 0)
%         err = [err, (res(p(1), p(2)) - tmp)];
% %         %res(~mask) = NaN;
%         hh = surf(res, 'EdgeColor', 'none');
% %         title(['Minimal surface - Jacobi method, ' int2str(iter) ...
% %             ' iterations, rastr [x, y] = [' int2str(rx) ', ' int2str(ry) '] with time: '...
% %             num2str(toc) 's, error in [0, 0] = ' num2str((res(pnum) - tmp))], 'Fontresize', 15);
%         hh.XData = linspace(min(lx), max(lx), rx);
%         hh.YData = linspace(min(ly), max(ly), ry);
%         view(-70, 40); drawnow;
%     end
end

endtime = toc

% figure; plot(linspace(0, iter, numel(err(2:end))), err(2:end), 'LineWidth', 3);
% % set(gca, 'XTickLabel', 0:1e3:iter);
% title('Delta progress function', 'FontSize', 16);
% xlabel('Iterations [-]', 'FontSize', 14); ylabel('Diference [-]', 'FontSize', 14);
% grid on

res = reshape(res, ry, rx);
res(~mask) = NaN;
figure; hh = surf(rot90(rot90(res)), 'EdgeColor', 'none'); %view(0, 90);
hh.XData = linspace(min(lx), max(lx), rx);
hh.YData = linspace(min(ly), max(ly), ry);

axmin = min(min(lx), min(ly));
axmax = max(max(lx), max(ly));
axis([axmin, axmax, axmin, axmax, min(min(res)), max(max(res))]);

hold on
plot3(max(lx)*lx/(scalepar), max(ly)*ly/(scalepar), lz, 'LineWidth', 5, 'Color', 'red');
hold off
title(['Convolution method, time: ' num2str(endtime*1e3) ' ms, iterations: ' num2str(iter)], 'FontSize', 15);
xlabel('x'); ylabel('y'); zlabel('z')





