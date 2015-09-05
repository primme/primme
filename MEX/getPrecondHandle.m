function y = getPrecondHandle(x)
%  PRIMME_MEX calls this function to perform the preconditioning operation.
%  If P1, P2 are matrices, it performs P2\(P1\x)) with the block vector x.
%  Otherwise, it passes the block vector x to the user defined P1(x) function 
%  handle to perform the preconditioning. The result y is returned to PRIMME_MEX.

    global P1;        % P1 is the first preconditioner matrix or function
    global P1matrix;  
    global P2;        % P2 is the second preconditioner matrix
    global eigsFunCallFlag; % mark if primme_egis is called by users
        
    % check user calls primme_eigs
    if  ~isempty(eigsFunCallFlag)
        if P1matrix % P1 is a matrix for primme_eigs
            y = P1\x;
            if  ~isempty(P2) % P2 exists and is a matrix for primme_eigs
                y = P2\y;
            end
        else % P1 is a matrix function for primme_eigs
            y = P1(x);
        end
    end
end

