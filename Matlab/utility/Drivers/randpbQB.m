% QB block algorithm, see http://arxiv.org/abs/1503.07157
function [Q,B] = randpbQB(A,q,s,kstep,nstep)

stepmult = nstep*kstep;
m   = size(A,1);
n   = size(A,2);
Q   = zeros(m,stepmult);
B   = zeros(stepmult,n);

Omega   = randn(n,stepmult);

for istep = 1:nstep
  ind = kstep*(istep-1) + (1:kstep);
  fprintf('istep =%d, ind(1) = %d, ind(end) = %d\n', istep, ind(1), ind(end));

  tstart = tic;
  Y   = A*Omega(:,ind);
  telapsed = toc(tstart);
  fprintf('elapsed time for building Y: %f sec\n', telapsed);


  tstart = tic;
  for j = 1:q
    if mod(2*j-2,s) == 0
        [Qnew,~] = qr(Y,0);
    else
        Qnew = Y;
    end
    Y          = A'*Qnew;
        
    if mod(2*j-1,s) == 0
        [Qnew,~] = qr(Y,0);
    else
        Qnew = Y;
    end
    Y          = A*Qnew;
  end
  [Qnew,~] = qr(Y,0);
  telapsed = toc(tstart);
  fprintf('elapsed time for power method: %f sec\n', telapsed);

  tstart = tic;
  if (istep > 1)
    J          = 1:((istep-1)*kstep);
    %fprintf('J(1) = %d, J(end) = %d\n', J(1), J(end));
    Y          = Qnew - Q(:,J)*(Q(:,J)'*Qnew);
    [Qnew,~] = qr(Y,0);
  end
  telapsed = toc(tstart);
  fprintf('elapsed time for reorthogonalization: %f sec\n', telapsed);

 
  tstart = tic;   
  Bnew     = Qnew'*A;
  telapsed = toc(tstart);
  fprintf('elapsed time for Bnew: %f sec\n', telapsed);

  tstart = tic;   
  A        = A - Qnew*Bnew;
  telapsed = toc(tstart);
  fprintf('elapsed time for update A: %f sec\n', telapsed);

  tstart = tic;   
  Q(:,ind) = Qnew;
  B(ind,:) = Bnew;
  telapsed = toc(tstart);
  fprintf('elapsed time for updating Q and B: %f sec\n', telapsed);

end


