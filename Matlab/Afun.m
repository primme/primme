function y = Afun(x,tflag)
 %   This example shows how to specify user' preconditioner function P and
 %   P'for primme_svds function. P may be a function handle PFUN such that 
 %   PFUN(X,'notransp') returns P\X and PFUN(X,'transp') returns P'\x.
    global primmeA;
    if strcmp(tflag,'notransp')
        y = primmeA*x;
    else
        y = primmeA'*x;
    end
end
