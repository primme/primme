function plotProjectionErrors(projerr)
%
% PROJERRPLOT(DATA)
%   DATA (i,j,k,l) struct
%       i - columns d
%       j - rank m
%       k - spectral decay zeta
%       l - sketch sizes ell
%       Fields:
%           1. left: left projection errors for us from FFD
%           2. right: right projection errors for vs from FFD
%     us,vs have ell > m columns. We consider only the m largest of them for our
%     projection error. It is possible that for k > m, the sketch projection
%     error is less than the exact projection error for m.
%
    close all;
    % Data dimensions
    w = size(projerr,1);   % no. of cols
    x = size(projerr,2);   % no. of ranks
    y = size(projerr,3);   % no. of noise
    z = size(projerr,4);   % no. of sketches
    n = projerr(1).rows;   % row sizes (default n = 100,000)
    d = projerr(1).cols;   % columns sizes
    m = projerr(1).ranks;
    ells = projerr(1).ells;
    z = numel(ells);
    
    % Get axes

    for i = 1:w
        for j = 1:x
            for k = 1:z
                l1(k) = projerr(i,j,2,k).left;
                l2(k) = projerr(i,j,3,k).left;
                r1(k) = projerr(i,j,2,k).right;
                r2(k) = projerr(i,j,3,k).right;
            end
            g1(i,j) = max(l1);
            g2(i,j) = max(l2);
            g3(i,j) = max(r1);
            g4(i,j) = max(r2);
            gl(i,j) = max(g1(i,j),g2(i,j));
            gr(i,j) = max(g3(i,j),g4(i,j));
        end
        Gl(i) = max(gl(:,i));
        Gr(i) = max(gr(:,i));
    end
%     G = max(Gl,Gr);
    

    
    % Plot labels
    for i = 1:numel(ells)
        axtix(i) = ells(i);
        axval{i} = i*10;
    end

    % Plot titles
    for i = 1:w
        for j = 1:x
%             top = G(j); uncomment if using 1 axis for left & right across
%             ranks
            for k = 1:y
%                 
                for p = 1:z
                    lpe(p) = projerr(i,j,k,p).left;
                    rpe(p) = projerr(i,j,k,p).right;
                end
                
                str = '$N$ = %d, $d$ = %d, $m$ = %d.';
                str = sprintf(str,n,d(i),m(j));
                
                % left projerr
                figure(10*i + (2*j-1));
                lplot = semilogy(ells,lpe,'LineWidth',2.5);
                hold on;
                xticks(axtix);
                xticklabels(axval);
                xlabel('Dimension of the FFD sketch','Interpreter','latex');
                ylabel('Projection Error','Interpreter','latex');
                leftplot(k) = lplot(1);
                lleg = legend(leftplot,{'$\zeta = 0.5$','$\zeta = (1+\sqrt{d/m})/2$','$\zeta = \sqrt{d/m}$','$\zeta = 10\sqrt{d/m}$'},'Interpreter','latex');
                title1 = {'$|| A - \tilde{U} \tilde{U}^T A ||_F / || A - U U^T A ||_F$';...
                    '$\tilde{U}$: $m$ approximate left singular vectors of the sketch for various spectral decays $\zeta$'};
                title1{end+1} = str;
                title(title1,'Interpreter','latex');
                top = Gl(j);
                margin = (top-1)*0.03;
                ylim([1-margin, top+margin]);
                

                % right projerr
                figure(10*i + 2*j);
                rplot = semilogy(ells,rpe,'LineWidth',2.5);
                hold on;
                xticks(axtix);
                xticklabels(axval);
                xlabel('Dimension of the FFD sketch','Interpreter','latex');
                ylabel('Projection Error','Interpreter','latex');
                rightplot(k) = rplot(1);
                rleg = legend(rightplot,{'$\zeta = 0.5$','$\zeta = (1+\sqrt{d/m})/2$','$\zeta = \sqrt{d/m}$','$\zeta = 10\sqrt{d/m}$'},'Interpreter','latex');
                title2 = {'$|| A - A \tilde{V} \tilde{V}^T ||_F / || A - A V V^T ||_F$';...
                    '$\tilde{V}$: $m$ approximate right singular vectors of the sketch for various spectral decays $\zeta$'};
                title2{end+1} = str;  
                title(title2,'Interpreter','latex');
                top = Gr(j);
                margin = (top-1)*0.03;
                ylim([1-margin, top+margin]);
            end
        end
    end
