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
