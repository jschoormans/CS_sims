function res = CoilMultiplicationit(a,b)

for nc=1:size(a.sens,3)
res(:,:,nc)=a.sens(:,:,nc).*b
end