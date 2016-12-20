function res = mtimes(a,b)

%-------------------------------------------------------------------------

if a.adjoint
    res = CoilMultiplicationit(a,b);
else
    res = CoilMultiplicationt(a,b);
end

