function  res = Wavelet_SQKP(size,dims)
%operator making square ksp with size of next 2^n

res.W = Wavelet('Daubechies',4,4);	% Wavelet
res.SQ=SQKSP(size,[1,2]); % SQUARE KSP OPERATOR
res.adjoint = 0;
res = class(res,'Wavelet_SQKSP');

