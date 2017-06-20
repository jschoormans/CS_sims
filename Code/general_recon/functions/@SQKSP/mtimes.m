function res = mtimes(a,b)


if a.adjoint
	res = sqkspadjoint(b,a.dims,a.origsize);

else
	res = sqkspforward(b,a.dims);

end




    
