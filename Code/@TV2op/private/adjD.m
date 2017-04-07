function res = adjD(y)

res = zeros(size(y,1),size(y,2));
res = adjDx(y(:,:,1)) + adjDy(y(:,:,2));

return;


function res = adjDy(x)
res = x(:,[1,1:end-1]) - 2.*x + x(:,[2:end,end]);
res(:,1) = -x(:,1);
res(:,end) = x(:,end-1);

function res = adjDx(x)
res = x([1,1:end-1],:) - 2.*x + x([2:end,end],:);
res(1,:) = -x(1,:);
res(end,:) = x(end-1,:);


