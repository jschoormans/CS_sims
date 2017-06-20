function res = CoilMultiplication(sens)

%-------------------------------------------------------------------------
% 2D coil multiplication operator (x,y,nc)
%-------------------------------------------------------------------------

res.adjoint = 0;
res.sens=sens;
res = class(res,'CoilMultiplication');
