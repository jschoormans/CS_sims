function [PSF,OTF]=PSFinF(N,xvals,amplitudes,visualize)
% 2D PSF/OTF for non-discrete values

assert(length(xvals)==length(amplitudes)); 

T=N;
omega0=2*pi/T;

k=[1:N]

E=zeros(1,128)

for i=1:length(xvals)
t0=xvals(i)
C=ones(1,N);
D=C.*exp(-1i*omega0.*k.*t0).*amplitudes(i)
E=E+D
end
OTF=(E);
OTF=ifftshift(E);
PSF=(ifft(OTF));
if visualize ==1
    
    f=linspace(-63,64,128)
    
    figure(1);
    subplot(221);
    stem(f,real(OTF)); title('real OTF')
    subplot(223)
    stem(f, imag(OTF)); title('imaginary OTF')
    subplot(222)
    plot([1:N],abs(PSF),'.-')  ; title('PSF') 
end

end
