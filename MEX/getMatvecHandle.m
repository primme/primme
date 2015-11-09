function y = getMatvecHandle(x)
%  PRIMME_MEX calls this function to perform the MATVEC. If primmeA is a matrix,
%  the MATVEC is performed by multiplication. Otherwise, it passes the block 
%  vector x to the user defined primmeA function handle to perform the MATVEC.
%  The result y is returned to PRIMME_MEX.

    global primmeA; % primmeA is a matrix or a matrix function
    global Amatrix; % mark if A is  a matrix or a matrix function
    global eigsFunCallFlag; % mark if primme_eigs is called by users
    
    % check user calls primme_eigs 
    if  ~isempty(eigsFunCallFlag) 
        if Amatrix % primmeA is a matrix for primme_eigs
            y = primmeA * x;
        else % primmeA is a matrix function for primme_eigs
            y = primmeA(x);
        end
    end
end

