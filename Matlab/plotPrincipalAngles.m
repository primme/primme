function plotPrincipalAngles(data)
%
% ANGLEPLOT(DATA,ELL)
%   DATA (i,j,k,l) struct
%       i - columns d
%       j - rank m
%       k - spectral decay zeta
%       l - sketch sizes ell
%       Fields:
%           1. Principal angles between the m sing.vectors of A and the ell-dim
%              space of the sketch
%               a. left
%               b. right
%
%   ELL -- ell = 10:10:200 -- default ell = 50
% 
%   REQUIRES violinplot.m & Violin.m
%
    close all;

    ell = 50;
    l = ell/10;
    w = size(data,1);   % no. of cols
    x = size(data,2);   % no. of ranks
    y = size(data,3);   % no. of noise
    z = size(data,4);   % no. of sketches
    n = data(1).rows;   % row sizes (default n = 100,000)
    d = data(1).cols;   % columns sizes
    m = data(1).ranks;

    colors = get(gca,'ColorOrder');
    mlabels = {'$m = 10$','$m = 20$','$m = 50$'};
    zlabels = {'0.5','1.724','2.450','24.495',...
        '0.5','1.724','2.450','24.495',...
        '0.5','1.724','2.450','24.495'};
    L = zeros(50,12); % max rank size 50; 3 four-column sets
    % R = zeros(50,12); % max rank size 50; 3 four-column sets

    for i = 1:w
        c = 1;
        for j = 1:x
            for k = 1:y
            % plot a fixed sketch size ell
                % since datasets have different # rows, violinplots.m requires 1
                % unified dataset so plots don't overlap

                L(1:m(j),c) = [data(i,j,k,l).left]; %left SVs full
                % R(1:m(j),c) = [data(i,j,k,l).right]; %right SVs full
                c = c + 1;
            end
        end
        
        str = '$N$ = %d, $d$ = %d, $ell$ = %d.';
        str = sprintf(str,n,d(i),ell);

        % plot pretty scatter plots using violinplots.m

        figure(2*i-1);
            hold on;
            vL = violinplot(L);
            VL = [vL(1).ScatterPlot vL(5).ScatterPlot vL(9).ScatterPlot];
            legvL = legend(VL,mlabels,'Interpreter','latex');
                legvL.FontSize = 10;
                legvL.LineWidth = .75;
                legvL.Location = 'southeast';
            set(gca,'Yscale','log');
            xticklabels(zlabels);
            xlabel('Spectral decay $\zeta$','Interpreter','latex');
            ylabel('Sines of the $m$ principal angles','Interpreter','latex');
            title1 = {'Subspace angles between the $m$ largest left singular vectors of';... 
                '$A = M + N/\zeta$ and the $ell$-dimensional left space of the sketch $B$.';...
                '$A$ $(N \times d)$, $B$ $(N \times m)$, $U_k$ $(N \times k)$, rank($M$) = $m$.'};
            title1{end+1} = str;
            title(title1,'Interpreter','latex');
            ax.XGrid = 'on';

        % figure(2*i);
        %     hold on;
        %     vR = violinplot(R);
        %     VR = [vR(1).ScatterPlot vR(5).ScatterPlot vR(9).ScatterPlot];
        %     legvR = legend(VR,mlabels,'Interpreter','latex');
        %         legvR.FontSize = 10;
        %         legvR.LineWidth = .75;
        %         legvR.Location = 'southeast';
        %     set(gca,'Yscale','log');
        %     xticklabels(zlabels);
        %     xlabel('Spectral decay $\zeta$','Interpreter','latex');
        %     ylabel('Sines of the $m$ principal angles','Interpreter','latex');
        %     title2 = {'Subspace angles between the $m$ largest right singular vectors of';... 
        %               '$A = M + N/\zeta$ and the $ell$-dimensional right space of the sketch $B$.';...
        %               '$A$ $(N \times d)$, $B$ $(N \times m)$, $U_k$ $(N \times k)$, rank($M$) = $m$.'};
        %     title2{end+1} = str;
        %     title(title2,'Interpreter','latex');        
        %     ax.XGrid = 'on';

        for a = 1:size(L,2)

            vL(a).WhiskerPlot = 'off';
            vL(a).EdgeColor = 'none';
            vL(a).ViolinPlot.FaceAlpha = 0;
            % vR(a).WhiskerPlot = 'off';
            % vR(a).EdgeColor = 'none';
            % vR(a).ViolinPlot.FaceAlpha = 0;

            % set colors for all 3 ranks
            if (a <= 4)
                vL(a).ViolinColor = colors(1,:);
                % vR(a).ViolinColor = colors(1,:);
            elseif (4 < a) && (a <= 8)
                vL(a).ViolinColor = colors(2,:);
                % vR(a).ViolinColor = colors(2,:);
            elseif (8 < a) && (a <= 12)
                vL(a).ViolinColor = colors(4,:);
                % vR(a).ViolinColor = colors(4,:);
            end
        end
        clear L;
        % clear R;
    end
end
