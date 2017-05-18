function Ki=makeUndersampledK(K,MNSA,sig)

Noise=randn(size(K)).*sig; %noise contribution to k-space;
Ki=K+Noise;
Ki(MNSA==0)=0; %remove zeroes;

end