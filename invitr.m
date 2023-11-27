function [x,iter] = invitr(A, ep,  numitr)
%INVITR Inverse iteration
%[x,iter] = invitr(A, ep,  numitr) computes an approximation x, smallest
%eigenvector using inverse iteration.  initial approximation is vector of ones,
%ep is the tolerance and numitr is the maximum number of iterations.
%If the iteration converged, iter is the number of iterations
%needed to converge.  If the iteration did not converge,
%iter contains numitr.
%This program implements Algorithm in
%http://www.netlib.org/utk/people/JackDongarra/etemplates/node96.html
%input  : Matrix A, ep and integer numitr
%output : vector x and integer iter

    [m,n] = size(A);
    if m~=n
        disp('matrix A  is not square');
        return;
    end;

    y=ones(n,1);

    for k = 1 :  numitr
        iter = k;
        v = y/norm(y,2);
        y = A\v;
        th =v'*y;
        if norm(y-th.*v,2) < ep*abs(th)
            break;
        end;
    end;
      x = y/th;
end

