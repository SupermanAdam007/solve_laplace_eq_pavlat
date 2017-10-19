function [b, maskbound, mask, scalepar] = getBound(lx, ly, lz, rx, ry)
%Function getBound makes from paratetric curve lx, ly, lz boundary condition
%for Laplace equation.
%
%   OUTPUT:
%   b...........projection 3D parametric curve to 2D plane (matrix of 
%               doubles - value means high in z axe)
%   maskbound...logical mask of places where is boundary condition
%   mask........logical mask: true-inside of bound, false-outside of bound
%   scalepar....scalar which show nessesery scale of curve for defined
%               rastr
%
%   INPUT:
%   lx, ly, lz..vectors of points which define parametric curve
%   rx, ry......size of rastr in x, y axe


b = zeros(ry, rx);
cx = round(size(b, 2)/2); cy = round(size(b, 1)/2);
scalepar = max(max(ly, lx));
% scalepar = 1;
lyr = round(ly*(cy-1)/scalepar); lxr = round(lx*(cx-1)/scalepar);
% maskbound = zeros(ry, rx);
mask = zeros(size(b));

for ii = 1:numel(lx)
%     coors(:,ii) = [cy - lyr(ii); cx - lxr(ii)];
    b(min(cy - lyr(ii), ry), min(cx - lxr(ii), ry)) = lz(ii);
end
maskbound = abs(b) > 0;

for ii = 1:size(b, 1)
    fromleft = 1; % jdu z leva
    while ((b(ii, fromleft) == 0) && (fromleft ~= size(b, 2))) % pocitam dokud nenarazim na nenulovej
        fromleft = fromleft + 1;
    end
    fromright = size(b, 2);
    while ((b(ii, fromright) == 0) && (fromleft ~= size(b, 2)))
        fromright = fromright - 1;
    end
    if (size(fromleft:fromright, 2) > 1)
        mask(ii, fromleft:fromright) = 1;
    end
end

for ii = 1:size(b, 2)
    fromtop = 1; % jdu z leva
    while ((b(fromtop, ii) == 0) && (fromtop ~= size(b, 1))) % pocitam dokud nenarazim na nenulovej
        fromtop = fromtop + 1;
    end
    fromdown = size(b, 1);
    while ((b(fromdown, ii) == 0) && (fromtop ~= size(b, 1)))
        fromdown = fromdown - 1;
    end
    if(size(1:fromtop, 2) > 1)
        mask(1:fromtop-1, ii) = 0;
    end
    if(size(size(b, 1):-1:fromdown, 2) > 1)
        mask(max(size(b, 1):-1:fromdown + 1, 1), ii) = 0;
    end
end
mask = logical(mask);

end
