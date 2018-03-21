function plotProjectionErrors(data)
%
% PROJERRPLOT(DATA)
%   DATA (i,j,k,l) struct
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
    close all;
    
    % dimensions of the struct
    w = size(data(1).rows,2);       % no. of cols
    x = size(data(1).ranks,1);      % no. of ranks
    y = size(data(1).zetas,2);      % no. of noise
    z = size(data(1).ells,2);       % no. of basis sizes
    
    % dimensions of the problem
    nn = data(1).rows;       % row sizes
    dd = data(1).cols;       % columns sizes
    mm = data(1).ranks;      % ranks
    zetas = data(1).zetas;  % noise 
    ells = data(1).ells;    % basis sizes
    
    % Get axes

    for i = 1:w
        for j = 1:x
            for k = 1:z
                if (data(i,j,2,k).ProjErr.left > data(i,j,3,k).ProjErr.left)
                    l(k) = data(i,j,2,k).ProjErr.left;
                else
                    l(k) = data(i,j,2,k).ProjErr.right;
                end
                if (data(i,j,2,k).ProjErr.right > data(i,j,3,k).ProjErr.right)
                    r(k) = data(i,j,2,k).ProjErr.right;
                else
                    r(k) = data(i,j,3,k).ProjErr.right;
                end
            end
            tmp = max(l);
            gl(i,j) = tmp;
            tmp = max(r);
            gr(i,j) = tmp;
        end
    end
    
    for j = 1:x
        Gl(j) = max(gl(:,j));
        Gr(j) = max(gr(:,j));
    end
    clear l r gl gr;

    
    % Plot labels
    for i = 1:numel(ells)
%         axtix(i) = ells(i);
        axval{i} = i*10;
    end

    % Plot
    for i = 1:w
    n = nn(i);
    d = dd(i);
            
        for j = 1:x
        m = mm(j,i);
%             top = G(j); uncomment if using 1 axis for left & right across
%             ranks
            for k = 1:y
%                 
                for p = 1:z
                    lpe(p) = data(i,j,k,p).ProjErr.left;
                    rpe(p) = data(i,j,k,p).ProjErr.right;
                end
                
                subtitle = '$N$ = %d, $d$ = %d, $m$ = %d.';
                subtitle = sprintf(subtitle,n,d,m);
                
                % left projerr
                figure(10*i + (2*j-1));
                lplot = semilogy(ells,lpe,'LineWidth',2.5);
                hold on;
                xticks(ells);
                xticklabels(axval);
                xlabel('Dimension of the FFD sketch','Interpreter','latex');
                ylabel('Projection Error','Interpreter','latex');
                leftplot(k) = lplot(1);
                lleg = legend(leftplot,{'$\zeta = 0.5$','$\zeta = (1+\sqrt{d/m})/2$','$\zeta = \sqrt{d/m}$','$\zeta = 10\sqrt{d/m}$'},'Interpreter','latex');
                title1 = {'$|| A - \tilde{U} \tilde{U}^T A ||_F / || A - U U^T A ||_F$';...
                    '$\tilde{U}$: $m$ approximate left singular vectors of the sketch for various spectral decays $\zeta$'};
                title1{end+1} = subtitle;
                title(title1,'Interpreter','latex');
                top = Gl(j);
                margin = (top-1)*0.03;
                ylim([1-margin, top+margin]);
                

                % right projerr
                figure(10*i + 2*j);
                rplot = semilogy(ells,rpe,'LineWidth',2.5);
                hold on;
                xticks(ells);
                xticklabels(axval);
                xlabel('Dimension of the FFD sketch','Interpreter','latex');
                ylabel('Projection Error','Interpreter','latex');
                rightplot(k) = rplot(1);
                rleg = legend(rightplot,{'$\zeta = 0.5$','$\zeta = (1+\sqrt{d/m})/2$','$\zeta = \sqrt{d/m}$','$\zeta = 10\sqrt{d/m}$'},'Interpreter','latex');
                title2 = {'$|| A - A \tilde{V} \tilde{V}^T ||_F / || A - A V V^T ||_F$';...
                    '$\tilde{V}$: $m$ approximate right singular vectors of the sketch for various spectral decays $\zeta$'};
                title2{end+1} = subtitle;  
                title(title2,'Interpreter','latex');
                top = Gr(j);
                margin = (top-1)*0.03;
                ylim([1-margin, top+margin]);
            end
        end
    end
