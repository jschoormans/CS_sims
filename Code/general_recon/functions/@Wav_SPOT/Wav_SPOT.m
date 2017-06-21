function  res = Wav_SPOT(size)

res.W = opWavelet2(size(1),size(2),'Haar');	% Wavelet
res.size=size; % SQUARE KSP OPERATOR
res.numel=size(1)*size(2);
res.adjoint = 0;
res = class(res,'Wav_SPOT');

