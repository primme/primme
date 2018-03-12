function [B,nSVD] = fastFD(A,k)
    n = size(A,1); d = size(A,2);
    k = 2*k;                % TEMP SKETCH SIZE
    B = zeros(k,d);
    indB = find(~any(B,2)); % INDEX ALL NON-ZERO ROWS OF B
    i = 1;                  % KEEP TRACK OF DATA SAMPLES APPENDED
    nSVD = 0;
    %% FAST FREQUENT DIRECTIONS ALGORITHM
    while i <= n
        % APPEND DATA
        if ~isempty(indB)
            % INSERT NEXT DATA SAMPLE INTO FIRST NON-ZERO ROW OF B
            B(indB(1),:) = A(i,:);
            indB(1) = [];
            i = i + 1;
        end
        % UPDATE SKETCH
        if isempty(indB)
            [~,S,V] = svd(B,'econ');
            nSVD = nSVD + 1;
            Sprime = reduceRankFast(S,k,1);
            B = Sprime*V';

            % INDEX REMAINING ALL-ZERO ROWS OF B
            indB = find(~any(B,2));
        end
    end
B = B(1:k/2,:);   % RETURN APPROX B OF SIZE k/2
end
