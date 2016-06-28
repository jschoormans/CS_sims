function Ku_N=addNoise(Ku,sigma)
Nr=randn(size(Ku)).*sigma;   %add real noise
Nim=randn(size(Ku)).*sigma*1i;   %add real noise


Ku_N=Ku+Nr+Nim;
M=Ku_N~=0;
Ku_N=Ku_N.*M;
end