function [kspace_sq] = squareksp(ksp);

S=size(ksp);
kspace_sq=ksp;
for ii=1:min(length(S),3); %only first 3 dimensions (spatial)
    d=zeros(length(S),1);
    if (log2(S(ii)))==ceil((log2(S(ii))));
    else
        d(ii)=(2^(ceil(log2(S(ii))))-S(ii))/2;
        if mod(d(ii),2)==0
              kspace_sq=padarray(kspace_sq,d,'both');
        else %pad one more on the left
                kspace_sq=padarray(kspace_sq,ceil(d),'pre');
                kspace_sq=padarray(kspace_sq,floor(d),'post');
        end
    end
end
end