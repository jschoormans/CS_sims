function [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra(Gx,Nx,Lx)
%inputs:
%Gx: gradient strength (in T) 
%Nx: number of pixels in FOV
%Lx: size of FOV (in m) 

%outputs: 
% PFCE,PFCE_alpha,PFOB,PFOB_alpha:locations (pixels) and intensities (A.U.) of chemical shift peaks  

%% FOR PFCE
PFCE=0;
PFCE_alpha=20;
%% FOR PFOB
df=[26.7,8.6,-27.1,-31.6,-36.2] %ppm

gamma=40.052e6  %gyromagnetic ratio for 19F [Hz/T]

dx=2*pi*df/(gamma*Gx)

PFOB=Nx*dx/Lx % in pixels
PFOB_alpha=[2,3,2,6,2]

figure; 
hold on
stem(PFCE,PFCE_alpha)
stem(PFOB,PFOB_alpha)
hold off
xlabel('pixels')
title('PSF of PFCD and PFOB')
end