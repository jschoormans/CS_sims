function [width] = SigmoidFitting(Y,Center,visualize,normalizeoption)
if ~exist('normalizeoption')
    normalizeoption=0; 
end
% Fits a Sigmoid function to an Image Edge 
 
X=[1:length(Y)];


% normalize Y(0 to 1?)
Yn=Normalize(Y,normalizeoption);

%%FIND CENTER (CAN ONLY BE DONE WITH ENOUGH SNR!)
if isempty(Center)
Center=FindCenter(Yn)
end
% LS FIT OF SIGMOID FUNCTION TO EDGE
a=nlsqSigmoid(Yn.',X,Center,visualize);
%%RISE LENGTH CALCULATION
width=abs(4.4/a);

end


function Center=FindCenter(Yn)

[ind1, val1]=find((Yn-0.5)>0,1,'first');

% if ind1==1
%     disp(['ind1 is 1!in1=',num2str(ind1),' Yn:',num2str(Yn)])
%     %remove first index (NOT SURE ABOUT THIS YET)
%     [val1,ind1]=find((Yn(2:end)-0.5)>0,1,'first');
%     ind1=ind1+1;
%     
%     if ind1==2
%         disp(['ind1 is 2!in1=',num2str(ind1),' Yn:',num2str(Yn)])
%     end
% end

% ind1=ind1+Center-6;
ind0=ind1-1;

val1=Yn(ind1)-0.5;
val0=Yn(ind0)-0.5;

Center=ind0+abs(val0)/(val1-val0);
end

function a=nlsqSigmoid(Yn,X,C,visualize)
a0=1; %starting point of nlsq
b0=1;
c0=0;
d0=C; %allow center to shift n pixels

lb=[-inf,-inf,-inf,C-5];
ub=[inf,inf,inf,C+5];
x0=[a0,b0,c0,d0]
fun= @(x)x(2)*(((1./(1+exp(-x(1).*(X-C))))))-Yn-x(3);
% fun= @(x)x(2)*(((1./(1+exp(-x(1).*(X-x(4)))))))-Yn-x(3);

x = lsqnonlin(fun,x0,lb,ub);
a=x(1);
if exist('visualize')
    if visualize==1
        figure(1);
        hold on
        plot(X-x(4),Yn,'r*')
        plot(X-x(4),fun(x)+Yn)
        hold off
        title('sigmoid fit')
    end
end
end

function Yn=Normalize(Y,normalizeoption)

if normalizeoption==0
% NORMALIZE BETWEEN 1 AND 0
Yn=(Y-min(Y))/(max(Y)-min(Y));

elseif normalizeoption==1
%NORMALIZE BETWEEN MEAN OF FIRST AND LAST 5 DATAPOINTS
Yn=(Y-mean(Y(1:5)))/(mean(Y(end-5:end))-mean(Y(1:5)));
else normalizeoption==2
    Yn=Y;
end


end


