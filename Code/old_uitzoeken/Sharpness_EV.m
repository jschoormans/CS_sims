function [Me]=Sharpness_EV(Image)
k=6;

N=length(Image); %what if Image is not square?


%normalize image by energy
I=Image./sqrt(sum(sum((Image).^2,1),2));
%computation of mean operator
mu=mean(I(:));
%computation of mean devaiation covariance operator
g=I-mu;
S=(1/(N*N-1))*(g*g.');
%computation of singular value matrix
[U,Sig,V]=svd(S);
%computation of eigenvalue matrix
Lambda=Sig.'*Sig;

%computation of proposed eigenvalue-based metric
Me=trace(Lambda(1:k,1:k));

end