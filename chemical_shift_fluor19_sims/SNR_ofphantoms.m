ii=(R'*K(1:length(K1)));

iii=abs(reshape(ii,[64 64])); 


noisepiece=iii(1:10,50:60);
imagepiece=iii(40:45,45:50);

figure(3); imshow(iii,[])
figure(1); imshow(noisepiece,[0 100]); 
figure(2); imshow(imagepiece,[0 100]);
noise=std(noisepiece(:))
signal=mean(imagepiece(:))

signal/noise