function res = mtimes(a,b)


if a.adjoint
	res = a.SQ'*(a.W'*b);

else
	res = a.W*(a.SQ*b);

end




    
