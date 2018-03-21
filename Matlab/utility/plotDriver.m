    function plotDriver(primme,svds,rsvd,fd)
%
%   PLOTDRIVER(PRIMME,SVDS,RSVD,FD)
%       Input:
%           PRIMME, SVDS, RSVD, FD data structures -- all info relevant to
%           parsing data and plotting is contained in fields
%               rows, cols, ranks, time, residuals, zetas, ells
%       Output:
%           Three plots comparing (in order) PRIMME, SVDS, RSVD, FD:
%           1. Runtime
%           2. Projection Error
%           3. Residuals
%
%
%   DOES THE SAME THING AS PLOTDRIVER2 EXCEPT IT ONLY PLOTS EVERY OTHER BASIS SIZE
%
    close all;

    % use SVDS for problem info that is common across all methods
    nn = svds(1).rows;
    dd = svds(1).cols;
    mm = svds(1).ranks;  % 3 x n: problem ranks stored in columns
    zz = svds(1).zetas;  % problem zetas stored in rows
%    ells = svds(1).ells; % row vector
%    ells = ells(1:2:end);

    % indexing parameters
    w = numel(nn);
    x = size(mm,1);
    y = size(zz,2);

    % legend labels
    leglabs = {'PRIMME','SVDS','RSVD','FD'};

    % colors
    colors = get(gca,'ColorOrder');

%--------------------------------------------------------------------------
    % Generate all plots in a single nested loop
count = 1;
for i = 1:w
n = nn(i);
d = dd(i);

    for j = 1:x
    m = mm(j,i);

        ells = fd(i,j,1,1).ells; % For FD we have many basis sizes ell
        z = numel(ells); 

        % common axis
        xtix = {};     % Note all methods use twice the shown space for working. 
        xtix{1} = m;    %PRIMME
        xtix{2} = '';
        xtix{3} = m;	%SVDS
        xtix{4} = '';
        xtix{5} = m+20;	%RSVD
        xtix{6} = '';
        for jj = 1:z     % FD many ells
            xtix{end+1} = ells(jj);
        end
        xtix{end+1} = '';

        for k = 1:y
        zeta = zz(i*j,k);

            for l = 1:z
                % runtime
                t_fd(l) = fd(i,j,k,l).time;
                % projection error
                err_fd(l) = fd(i,j,k,l).ProjErr.left;
                % residuals
                tmp = min(m, ells(l));
                res_tmp = fd(i,j,k,l).residuals';
                res_fd(1:tmp,l) = res_tmp(1:tmp);
            end

	    % Runtime
            t_primme = primme(i,j,k).time; 
	    t_svds = svds(i,j,k).time;
	    t_rsvd = rsvd(i,j,k).time;
	    % Residuals
            res_primme = primme(i,j,k).residuals';
	    res_svds = svds(i,j,k).residuals;
	    res_rsvd = rsvd(i,j,k).residuals';
	    % Projection Error
            err_primme = primme(i,j,k).ProjErr.left;
            err_rsvd = rsvd(i,j,k).ProjErr.left;

            T = [t_primme   0 t_svds 0 t_rsvd   0 t_fd];
            E = [err_primme 0 1      0 err_rsvd 0 err_fd];
            tmp = zeros(tmp,1);
            size(tmp)
            size(res_primme)
            size(res_rsvd)
            size(res_svds)
            size(res_fd)
            R = [res_primme tmp res_svds tmp res_rsvd tmp res_fd];

            clear t_primme t_svds t_rsvd t_fd;
            clear err_rsvd err_fd;
            clear res_primme res_svds res_rsvd res_fd;
            clear tmp;

            % looping subtitle string
            subtitle = '$n = %d$, $d = %d$, $m = %d$, $\\zeta = %.3f$.';
            subtitle = sprintf(subtitle,n,d,m,zz(j,k));

            % aspect ratio
            arx = 1.8;
            ary = 1.0;

            % Plot runtime
            figure(count);
            [plot_time,ind_f] = dataBarPlot(T,z);
            plot_time(1).CData = colors(7,:);
            plot_time(3).CData = colors(4,:);
            plot_time(5).CData = colors(1,:);
            for p = 1:z
                plot_time(ind_f(p)).FaceColor = colors(2,:);
            end
            for p = 1:numel(T)
                plot_time(p).FaceAlpha = .75;
                plot_time(p).EdgeAlpha = 0;
                plot_time(p).BarWidth = 0.8;
            end
            set(gca,'XTick',[1:5+2*numel(ells)]);
            set(gca,'XTickLabels',xtix);
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'PlotBoxAspectRatio',[arx,ary,1]);
            ax = gca;
            ax.XAxis.FontSize = 9;
            xlabel('Basis size used','Interpreter','latex');
            ylabel('Runtime (s)','Interpreter','latex');
            LT = [plot_time(1) plot_time(3) plot_time(5) plot_time(ind_f(1))];
            leg_time = legend(LT,leglabs);
            leg_time.Location = 'southeast';
            headtitle = {'Runtime of each method for various basis sizes, $l$'};
            headtitle{end+1} = subtitle;
            title(headtitle,'Interpreter','latex');
            outFile = ['runtime_',num2str(n),'_',num2str(d),'_',num2str(m),...
                '_',num2str(k),'.pdf'];
            saveas(gcf,outFile);clf;

            count = count + 1;

            % Plot projection error
            figure(count);
            [plot_projerr,ind_f] = dataBarPlot(E,z);
            plot_projerr(1).CData = colors(7,:);
            plot_projerr(3).CData = colors(4,:);
            plot_projerr(5).CData = colors(1,:);
            for p = 1:z
                plot_projerr(ind_f(p)).CData = colors(2,:);
            end
            for p = 1:numel(E)
                plot_projerr(p).FaceAlpha = .75;
                plot_projerr(p).EdgeAlpha = 0;
                plot_projerr(p).BarWidth = 0.8;
            end
            top = max(E);
            margin = (top-1)*0.03;
            ylim([1-margin, top+margin]);
            set(gca,'XTick',[1:5+2*numel(ells)]);
            set(gca,'XTickLabels',xtix);
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'PlotBoxAspectRatio',[arx,ary,1]);
            ax = gca;
            ax.XAxis.FontSize = 9;
            xlabel('Basis size used','Interpreter','latex');
            ylabel('$||A-U_mU_m^TA||$','Interpreter','latex');
            LP = [plot_projerr(1) plot_projerr(3) plot_projerr(5) plot_projerr(ind_f(1))];
            leg_projerr = legend(LP,leglabs);
            leg_projerr.Location = 'northeast';
            headtitle = {'Left projection errors of each method for various basis sizes, $l$'};
            headtitle{end+1} = subtitle;
            title(headtitle,'Interpreter','latex');
            outFile = ['projerr_',num2str(n),'_',num2str(d),'_',num2str(m),...
                '_',num2str(k),'.pdf'];
            saveas(gcf,outFile);clf;

            count = count + 1;

            % plot residuals
            figure(count);
            [plot_res,ind_f] = dataDotPlot2(R,z);
            plot_res(1).CData = colors(7,:);
            plot_res(3).CData = colors(4,:);
            plot_res(5).CData = colors(1,:);
            for p = 1:z
                plot_res(ind_f(p)).CData = colors(2,:);
            end
            VL = [plot_res(1) plot_res(3) plot_res(5) plot_res(ind_f(1))];
            leg_res = legend(VL,leglabs,'Interpreter','latex');
            leg_res.Location = 'southeast';
            set(gca,'Yscale','log');
            set(gca,'XTick',[1:5+2*numel(ells)]);
            set(gca,'XTickLabels',xtix);
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'PlotBoxAspectRatio',[arx,ary,1]);
            ax = gca;
            ax.XAxis.FontSize = 9;
            xlabel('Basis size used','Interpreter','latex');
            ylabel('$||A^T \tilde{u}_i - \tilde{\sigma}_i \tilde{v}_i||$','Interpreter','latex');
            headtitle = {'Residual norms of the largest $m$ singular triplets returned by each method';...
            'for various basis sizes, $l$'};
            headtitle{end+1} = subtitle;
            title(headtitle,'Interpreter','latex');
            outFile = ['res_',num2str(n),'_',num2str(d),'_',num2str(m),...
                '_',num2str(k),'.pdf'];
            saveas(gcf,outFile);clf;

            close all;
        end
    end
end
    end

function [barplot,index_fd] = dataBarPlot(data,count)
%
%     PLOTRUNTIME(DATA,COUNT)
%
%       Plots bar graph comparing PRIMME, SVDS, RSVD, FD for a single
%       problem size across all basis sizes
%
%       COUNT: number of bases for FD
%
    index_fd = 7;
    index_fd = index_fd:index_fd + count - 1;
    barplot = [];
    for i = 1:numel(data)
        b = bar(i,data(i),'FaceColor','flat');
        hold on;
        barplot = [barplot b];
    end
end
function [dotplot,index_fd] = dataDotPlot2(data,count)
%
%     DATABARPLOT(DATA,COUNT)
%
%       Plots dot graph comparing PRIMME, SVDS, RSVD, FD for a single
%       problem size across all basis sizes
%
%       COUNT: number of bases for RSVD & FD
%
%     x = [1:size(data,2)]';
    y = size(data,2);
    sz = 90;
    index_fd = 7;
    index_fd = index_fd:index_fd + count - 1;
    for i = 1:y
        x = i*ones(size(data,1),1);
        dotplot(i) = scatter(x,data(:,i),sz,'filled','MarkerFaceAlpha',0.4);
        hold on;
    end
end
