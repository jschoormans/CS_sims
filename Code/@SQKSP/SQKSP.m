function  res = SQKP(size,dims)
%operator making square ksp with size of next 2^n

res.adjoint = 0;
res.dims=dims;
res.origsize=size;
res = class(res,'SQKSP');

