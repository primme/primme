function y = getMatvecHandle(x)
%  PRIMME_MEX calls this function to get corresponding MATVEC function 
%  handle, then passes the blocked vector x to right MATVEC function and 
%  return the result y of MATRIX-VECTOR operations to PRIMME_MEX.  
%  Detailed explanation goes here

    global primmeA; % primmeA is a matrix or a matrix function
    global Amatrix; % mark if A is  a matrix or a matrix function
    global eigsFunCallFlag; % mark if primme_egis is called by users
    
    % check user calls primme_eigs 
    if  ~isempty(eigsFunCallFlag) 
        if Amatrix % primmeA is a matrix for primme_eigs
            y = primmeA * x;
        else % primmeA is a matrix function for primme_eigs
            y = primmeA(x);
        end
    end
end

