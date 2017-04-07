function  res = TV2op()

%res = TVOP()
%
% Implements a spatial 2nd order finite-differencing operator.

res.adjoint = 0;
res = class(res,'TV2op');

