function [x] = conjgrad_fluor(A, b, x,niter,realI)
    r = b - A * x;
    p = r;
    rsold = r' * r;
    
    f=figure;
    
    for i = 1:niter
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        
        imshow(reshape(abs(x),[64 64]),[]); drawnow; 
        mse=sum((x-realI).^2); 
        pause(1); 
        
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;        
        disp(['iter:',num2str(i),' MSE: ',num2str(mse)])

    end
end
