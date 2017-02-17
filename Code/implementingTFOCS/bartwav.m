function y= bartwav(x,n1,n2,mode); 


if mode==1
    y=bart('cdf97 3', reshape(x,n1,n2));
    y=reshape(y,n1*n2,1);
elseif mode==2
    y=bart('cdf97 3 -i', reshape(x,n1,n2));
    y=reshape(y,n1*n2,1);
else
    error('mode should be either 1 or 2')
end