function [B,nSVD] = fastFD(A,k)
    n = size(A,1); d = size(A,2);
    k2 = 2*k;                % TEMP SKETCH SIZE is twice
    B = zeros(k2,d);
    B(1:k,:) = A(1:k,:);    % Copy the first chunk
    i = k+1;                  % KEEP TRACK OF DATA SAMPLES APPENDED
    nSVD = 0;
    %% FAST FREQUENT DIRECTIONS ALGORITHM
    while i <= n
        % APPEND k more rows 
	upto = min(n,i+k-1); chunksize = upto-i+1;
        B(k+1:k+chunksize,:) = A(i:upto,:);
	i = upto+1;

        % UPDATE SKETCH
        [~,S,V] = svd(B(1:k+chunksize,:),'econ');
        nSVD = nSVD + 1;
        Sprime = reduceRankFast(S,k2,1);
	% There is little difference between the three choices below. Pick the middle one:
	%B = Sprime*V';
        B(1:k,:) = Sprime(1:k,1:k)*V(:,1:k)';
	%for j=1:k
	%	B(j,:) = Sprime(j,j)*V(:,j)';
	%end
    end
B = B(1:k,:);   % RETURN APPROX B OF SIZE k
end
%% REDUCE RANK
function Sprime = reduceRankFast(S,k,alpha)
    s = diag(S);
    sprime = zeros(size(s));

    skip = floor(k*(1-alpha)) + 1;
    if skip > 1
        sprime(1:skip) = s(1:skip);
    end

    dirac_ind = k - floor(k*alpha/2) + 1;
    if (skip < k) && (dirac_ind <= k)
        dirac = s(dirac_ind)^2;
        sprime(skip:end) = sqrt( max(s(skip:end).^2 - dirac,0) );
    end

    Sprime = diag(sprime);
end

