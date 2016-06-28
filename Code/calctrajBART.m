function [kspace, coordsdyns]=calctrajBART(K);  
%function to transform a Cartesian k-space with multiple dynamics into one list of sampled k-lines
%undersampled lines will be removed

%input:K(fe,pe1,pe2,nc,ndyn)
%output: coordinates and reshaped k-space with only sampled points (BART
%toolbox input)

%3-6-2016: TO DO: for even numbered signals, the trajectory numbers do not
%yet add up!

for nc=1:size(K,4);
    kspacedyns=[];
    coordsdyns=[];
    for ndyn=1:size(K,5);
    [k, coords]= kspacereshape(K(:,:,:,nc,ndyn));
    kspacedyns=[kspacedyns,k];
    coordsdyns=cat(3,coordsdyns,coords);
    end
    kspace(:,:,:,nc)=kspacedyns;
end



function [kspace, coords]= kspacereshape(K);
    n1=size(K,1);
    n2=size(K,2);
    n3=size(K,3);
    
    %calculate relevant k-space points and reshape
    Kr=reshape(K,[n1,n2*n3]);           %matrix to vector
    [index]=find(sum(Kr,1)~=0);         %non-sampled points
    kspace=Kr(:,index);
    
    I2=ones(1,n3)'*linspace(-1,1,n2);   %y-direction
    I3=linspace(-1,1,n3)'*ones(1,n2);   %z-direction
    I=I2+1j*I3;
    Cr=I(:);
    
    traj=Cr(index);

    coords=zeros([3,size(kspace)]);
    coords(1,:,:)=repmat(linspace(-1,1,n1)'.*n1,[1 length(index)]);
    coords(2,:,:)=ones(1,n1)'*real(traj).'.*n2;
    coords(3,:,:)=ones(1,n1)'*imag(traj).'.*n3;
end
end
