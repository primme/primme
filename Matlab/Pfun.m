function y = Pfun(x,model)
 % This example shows how to specify user' preconditioner function directly
 % for AHA, AAH, and aug. We use RIF to generate a preconditioner for AHA.
    global LL;
    global primmeA;
    [m,n] = size(primmeA);
    if strcmp(model, 'AHA') % we have RIF for A'A = LL*LL'
        y = eye(n)\x; % y = LL^-T * LL^-1 * x
    elseif strcmp(model, 'AAH')  % we have RIF for AA' = LL*LL'
        y = LL'\(LL\x); % y = LL^-T * LL^-1 * x
    else
        y = eye(m+n)\x; % apply identity matrix on aug
    end
end

