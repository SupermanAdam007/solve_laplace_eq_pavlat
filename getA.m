function A = getA(rx, ry)
%Function getA generate matrix size rx
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


xpart = sparse(2:(rx+2)-2, 1:(rx+2)-3, 1/4, rx, rx); %sparse: squeezing out every zero element
Ax = xpart + xpart';

ypart = sparse(2:(ry+2)-2, 1:(ry+2)-3, 1/4, ry, ry);
Ay = ypart + ypart';

A = kron(Ay, speye(rx)) + kron(speye(ry), Ax); % speye := sparse(eye(.., ..))

end

