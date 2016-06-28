function S= GenSensMaps(nx,ny,nS) %generate 2d sense maps

S=zeros(nx,ny,1,nS);

for n=1:floor(nS/2)
    for x=1:nx
        for y=1:ny
            S(x,:,n)=((x)^2+(y)^2)/()
        end
    end
    
end 
for n=floor(nS/2)+1:nS
    
end