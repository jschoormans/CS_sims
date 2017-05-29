    %find PSF by blind deconvolution
    
    linear_recon=opInverse(F2)*data; %linear recon of data 

    Spectrum_e=mean(abs(linear_recon_reshape(:,129)),2).';
    Spectrum_e=Spectrum_e./sum(Spectrum_e(:));
    Spectrum_e=Spectrum_e.*(Spectrum_e>0.03); % empirical spectrum
    
    [J P] = deconvblind(abs(linear_recon_reshape(1:128,129:133)),Spectrum_e.',10)
    offset=-min(round(PFOB));
    [J2 P2] = deconvblind(abs(linear_recon_reshape(10:20,1:128)),Spectrum_e,10)
    figure(7);hold on; plot(P); plot(P2); stem(Spectrum); hold off
    Spectrum=circshift(P.',[0 -17]);

    