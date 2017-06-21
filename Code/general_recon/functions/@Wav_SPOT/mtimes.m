function res = mtimes(a,b)


if a.adjoint
	res = reshape(a.W'*b(:),a.size);

else
	res = a.W*b(:);

end




    
