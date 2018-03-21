function plotPrincipalAngles(data,ell,R)
%
% PLOTPRINCIPALANGLES(DATA,ELL)
%   DATA (i,j,k,l) -- struct
%       i - columns d
%       j - rank m
%       k - spectral decay zeta
%       l - sketch sizes ell
%       Fields:
%           1. rows
%           2. cols
%           3. ranks
%           4. zetas
%           5. time
%           6. sine
%               a. left
%               b. right
%           7. ProjErr
%   
%   ELL -- Optional. Scalar input, basis size. Defaults to ell = 10 if no 
%          input is provided or if that problem size was never computed.
% 
%   R -- Optional. Basis size ell *must* be specified. If 'R' is specified, 
%        sines of right space will be plotted. Sines of left space always 
%        plotted by default.
% 
%   REQUIRES violinplot.m & Violin.m
% 
    close all;
    
    DO_L = 1; % Always print sines of left space by default
    DO_R = 0; % Don't print sines of right space by default
    
    if (nargin == 3) && (R == 'R')
        DO_R = 1;
    end
    if (nargin < 2)
        ell = 10;
    end
    
    % dimensions of the struct
    w = size(data(1).rows,2);       % no. of cols
    x = size(data(1).ranks,1);      % no. of ranks
    y = size(data(1).zetas,2);      % no. of noise
    z = size(data(1).ells,2);       % no. of basis sizes
    
    % dimensions of the problem
    n = data(1).rows;       % row sizes
    d = data(1).cols;       % columns sizes
    m = data(1).ranks;      % ranks
    zetas = data(1).zetas;  % noise 
    ells = data(1).ells;    % basis sizes
    
    % other struct related specs
    if isempty(find(ells == ell, 1))
        fprintf('\nNo data found for input basis size. Default set to minimum basis size ell = %d.\n',ells(1));
        ell = ells(1);  % default min basis size
        l = 1;          % index of default basis size
    else
        l = find(ells == ell);          % get index for ell we wish to plot
        fprintf('\nGenerating figures of principal angles for basis size ell = %d.\n',ell);
    end
    
    colors = get(gca,'ColorOrder'); % needed for uniformity across plots
    dim = x*y;

    for i = 1:w
        
        % labels for x axis
        for j = 1:size(zetas,1)/size(n,2)
            zetas(j,:) = round(zetas(j,:),3);
            zlabels{j} = zetas(j,:);
        end
        
        % header subtitle
        subtitle = '$N$ = %d, $d$ = %d, $ell$ = %d.';
        subtitle = sprintf(subtitle,n(i),d(i),ell);
        
        % legend labels
        mlabels = {};
        for j = 1:size(m,1)
            labs = 'm = %d';
            labs = sprintf(labs,m(j,i));
            mlabels{end+1} = labs;
        end 
        
        c = 1;
        L = zeros(max(m(:,i)),dim);
        R = zeros(max(m(:,i)),dim);
        for j = 1:x
            for k = 1:y
                % since datasets have different # rows, violinplots.m requires 1
                % unified dataset so plots don't overlap
                
                sinesFound = numel(data(i,j,k,l).sine.left);
                sines2Plot = min(m(j,i),sinesFound);
                L(1:sines2Plot,c) = [data(i,j,k,l).sine.left]; %left SVs full
                
                sinesFound = numel(data(i,j,k,l).sine.right);
                sines2Plot = min(m(j,i),sinesFound);
                R(1:sines2Plot,c) = [data(i,j,k,l).sine.right]; %right SVs full
                
                c = c + 1;
            end
        end
        
        % Plot pretty scatter plots using violinplots.m -------------------
        % Plot sines of left space
        if (DO_L)
            figure(2*i-1);
            hold on;
            vL = violinplot(L);
            % tune the appearance of the plots
            VL = [vL(1).ScatterPlot vL(5).ScatterPlot vL(9).ScatterPlot];
            legvL = legend(VL,mlabels,'Interpreter','latex');
            legvL = legend(VL,mlabels);
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
            title1{end+1} = subtitle;
            title(title1,'Interpreter','latex');
            ax.XGrid = 'on';
            for a = 1:size(L,2)
                vL(a).WhiskerPlot = 'off';
                vL(a).EdgeColor = 'none';
                vL(a).ViolinPlot.FaceAlpha = 0;

                % set colors for all 3 ranks
                if (a <= 4)
                    vL(a).ViolinColor = colors(1,:);
                elseif (4 < a) && (a <= 8)
                    vL(a).ViolinColor = colors(2,:);
                elseif (8 < a) && (a <= 12)
                    vL(a).ViolinColor = colors(4,:);
                end
            end
            
        L = [];
        end
        
        % Plot sines of right space
        if (DO_R)
            figure(2*i);
            hold on;
            vR = violinplot(R);
            VR = [vR(1).ScatterPlot vR(5).ScatterPlot vR(9).ScatterPlot];
            legvR = legend(VR,mlabels,'Interpreter','latex');
                legvR.FontSize = 10;
                legvR.LineWidth = .75;
                legvR.Location = 'southeast';
            set(gca,'Yscale','log');
            xticklabels(zlabels);
            xlabel('Spectral decay $\zeta$','Interpreter','latex');
            ylabel('Sines of the $m$ principal angles','Interpreter','latex');
            title2 = {'Subspace angles between the $m$ largest right singular vectors of';... 
                      '$A = M + N/\zeta$ and the $ell$-dimensional right space of the sketch $B$.';...
                      '$A$ $(N \times d)$, $B$ $(N \times m)$, $U_k$ $(N \times k)$, rank($M$) = $m$.'};
            title2{end+1} = subtitle;
            title(title2,'Interpreter','latex');        
            ax.XGrid = 'on';
            
            for a = 1:size(R,2)
                vR(a).WhiskerPlot = 'off';
                vR(a).EdgeColor = 'none';
                vR(a).ViolinPlot.FaceAlpha = 0;

                % set colors for all 3 ranks
                if (a <= 4)
                    vR(a).ViolinColor = colors(1,:);
                elseif (4 < a) && (a <= 8)
                    vR(a).ViolinColor = colors(2,:);
                elseif (8 < a) && (a <= 12)
                    vR(a).ViolinColor = colors(4,:);
                end
            end
            
        R = [];
        end        
    end
end
