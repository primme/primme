function y = getPrecondHandle(x)
%  PRIMME_MEX calls this function to get corresponding precondition function 
%  handle, then passes the blocked vector x to right preconditioning  
%  function and return the result y of MATRIX-VECTOR operations to PRIMME_MEX.  
%  Detailed explanation goes here

    global P1;        % P1 is the first preconditioner matrix or function
    global P1matrix;  
    %     P1matrix = 0, First preconditioner is a matrix funciton
    %     P1matrix = 1, First preconditioner is matrix for A
    %     P1matrix = 2, First preconditioner is matrix directly for ATA or OAAO                
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

