% [INPUT]
% data = A cell array of structures returned by the "execute_tests" function when the "tab" parameter is set to "false".
% det  = A boolean indicating whether to create a detailed plot for each time series of fund returns in the dataset (optional, default=false).

function plot_data(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'struct'},{'nonempty'}));
        p.addOptional('det',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    det = res.det;

    plot_data_internal(data,det);

end

function plot_data_internal(data,det)

    plot_market_data(data);
    plot_style_factors(data.DatesNum,data.SF);
    plot_funds_overview(data);
    
    if (det)
        bm = data.BM;
        bm_cum = cumprod(1 + bm) - 1;
        
        for i = 1:data.Frms
            plot_fund_details(data.DatesNum,bm,bm_cum,data.FrmsRet(:,i),data.FrmsNam{i},data.Grps(i));
        end
    end

end

function plot_fund_details(date,bm,bm_cum,ret,frm,grp)

    ret_ann = ((prod(1 + ret) ^ (12 / numel(ret))) - 1) * 100;
    ret_avg = mean(ret) * 100;
    ret_com = cumprod(1 + ret);
    ret_cum = ret_com - 1;
    ret_cum_max = cummax([1; ret_com]);
    ret_dd = (ret_com ./ ret_cum_max(2:end)) - 1;
    ret_dd_max = abs(min(ret_dd)) * 100;

    fig = figure();
    set(fig,'Name',[frm ' Details'],'Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 36]);
    plot(date,ret,'Color','b');
    hold on;
        plot(date,bm,'Color','r');
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor','b','ShowBaseLine','off');
	hold off;
    set(sub_1,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t1 = title(sub_1,'Returns','Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);
    
    sub_2 = subplot(13,9,[46 72]);
    plot(date,ret_cum,'Color','b');
    hold on;
        plot(date,bm_cum,'Color','r');
	hold off;
    set(sub_2,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t2 = title(sub_2,'Cumulative Returns','Units','normalized');
    t2_pos = get(t2,'Position');
    set(t2,'Position',[0.4783 t2_pos(2) t2_pos(3)]);
    
    sub_3 = subplot(13,9,[82 108]);
    plot(date,ret_dd,'Color','b');
    set(sub_3,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t3 = title(sub_3,sprintf('Drawdowns (Maximum: %.2f%%)',ret_dd_max),'Units','normalized');
    t3_pos = get(t3,'Position');
    set(t3,'Position',[0.4783 t3_pos(2) t3_pos(3)]);
    
    l = legend([l1 l2],'Benchmark',frm,'Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('%s (Group %d)\nAverage Return: %.2f%% | Annualized Return: %.2f%%',frm,grp,ret_avg,ret_ann));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_funds_overview(data)

    clrs = [
        0.984 0.502 0.447;
        0.553 0.827 0.780;
        1.000 0.882 0.435;
        0.745 0.729 0.855;
        0.502 0.694 0.827;
        0.992 0.706 0.384;
        0.702 0.871 0.412;
        0.988 0.804 0.898;
        0.851 0.851 0.851;
        0.737 0.502 0.741;
        0.800 0.922 0.773
	];

    f = data.Frms;
    n = data.Obs;
    
    bm = data.BM;
    bm_dn_idx = bm <= 0;
    bm_dn = sum(bm(bm_dn_idx));
    bm_up_idx = bm > 0;
    bm_up = sum(bm(bm_up_idx));

    grps = data.Grps;
    grps_uni = sort(unique(grps));
    grps_len = numel(grps_uni);
    grps_leg = gobjects(grps_len,1);

    cdn = [1; zeros(f,1)];
    cup = [1; zeros(f,1)];
    
    for i = 1:f
        ret = data.FrmsRet(:,i);

        cdn(i+1) = sum(ret(bm_dn_idx)) / bm_dn;
        cup(i+1) = sum(ret(bm_up_idx)) / bm_up;
    end

    cdn_x_max = (max(cdn) + 0.3);
    cdn_x_min = (min(cdn) - 0.3);
    cdn_y_max = (max(cup) + 0.3);
    cdn_y_min = (min(cup) - 0.3);
    
    cr_max = max(cdn_x_max,cdn_y_max);
    cr_min = min(cdn_x_min,cdn_y_min);
    
    ret_exc = data.FrmsRet - repmat(data.RF,1,f);
    ret_exc_ann = (prod(1 + ret_exc) .^ (12 / n)) - 1;
    ret_std_ann = sqrt(12) .* std(data.FrmsRet);
    sr = ret_exc_ann ./ ret_std_ann;
    sr_avg = mean(sr);

    fig = figure();
    set(fig,'Name','Funds Overview','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    b = boxplot([data.BM data.FrmsRet],'Symbol','k.');
    hold on;
        set(findobj(sub_1,'Type','Line','Tag','Box'),'Color','k');
        set(findobj(sub_1,'Type','Line','Tag','Median'),'Color','k');
        set(findobj(sub_1,'-regexp','Tag','\w*Whisker'),'LineStyle','-');
        h = get(b(5,:),{'XData','YData'});
        patch(h{1,1},h{1,2},clrs(1,:));       
        for i = 2:size(h,1),patch(h{i,1},h{i,2},clrs(grps(i-1)+1,:)),end
        chil = get(sub_1,'Children');
        set(sub_1,'Children',chil([end 1:end-1]));
        for i = 1:grps_len;grps_leg(i)=area(0,NaN,'FaceColor',clrs(i+1,:),'ShowBaseLine','off');end
    hold off;
    set(sub_1,'XTickLabel',[{'Benchmark'} data.FrmsNam]);
    title(sub_1,'Boxes');
    
    sub_2 = subplot(13,9,[73 112]);
    scatter(cdn,cup,36,[clrs(1,:); clrs(grps+1,:)],'filled');
    hold on;
        line([cr_min cr_max],[1 1],'Color',clrs(1,:));
        line([1 1],[cr_min cr_max],'Color',clrs(1,:));
        text(1.03,1.07,'Benchmark','Color',clrs(1,:),'Clipping','on');
        for i = 1:f;text(cdn(i+1)+0.03,cup(i+1)+0.03,data.FrmsNam{i},'Clipping','on');end
    hold off;
    grid on;
    set(sub_2,'Box','on');
    set(sub_2,'XLim',[cr_min cr_max]);
    set(sub_2,'YLim',[cr_min cr_max]);
    xlabel('Downside');
    ylabel('Upside');
    title(sub_2,'Capture Ratios');
    
    sub_3 = subplot(13,9,[78 117]);
    hold on;
        for i = 1:grps_len
            grp = grps_uni(i);
            sr_curr = sr;
            sr_curr(grps ~= grp) = NaN;
            bar(1:f,sr_curr,'FaceColor',clrs(grp+1,:));
        end
        text(1:f,sr,sprintfc('%.2f',sr),'HorizontalAlignment','center','VerticalAlignment','bottom');
        line([0.5 (f + 0.5)],[sr_avg sr_avg],'Color','k');
    hold off;
    set(sub_3,'XLim',[0.5 (f + 0.5)],'XTick',1:f,'XTickLabel',data.FrmsNam);
    title(sub_3,'Sharpe Ratios');

    if (grps_len > 1)
        l = legend(grps_leg,sprintfc('Group %d',grps_uni),'Location','southoutside','Orientation','horizontal');
        l_pos = get(l,'Position');
        set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.47 l_pos(3) l_pos(4)]);
    end

    figure_title(sprintf('Returns Overview\nObservations: %d | First: %s | Last: %s',data.Obs,datestr(data.DatesNum(1),'mm/yyyy'),datestr(data.DatesNum(end),'mm/yyyy')));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_market_data(data)

    fig = figure();
    set(fig,'Name','Market Data','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[19 58]);
    plot(data.DatesNum,data.BM,'Color','r');
    set(sub_1,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Returns');
	t1 = title(sub_1,sprintf('Benchmark\nAverage: %.2f%% | Annualized: %.2f%%',(data.BMAvg * 100),(data.BMAnn * 100)),'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[t1_pos(1) 1.1 t1_pos(3)]);
    
    sub_2 = subplot(13,9,[73 112]);
    plot(data.DatesNum,(cumprod(1 + data.BM) - 1),'Color','r');
    set(sub_2,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Cumulative Returns');

    sub_3 = subplot(13,9,[24 63]);
    plot(data.DatesNum,data.RF,'Color','b');
    set(sub_3,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Rates');
	t3 = title(sub_3,sprintf('Risk-Free Rate\nAverage: %.2f%% | Annualized: %.2f%%',(data.RFAvg * 100),(data.RFAnn * 100)),'Units','normalized');
    t3_pos = get(t3,'Position');
    set(t3,'Position',[t3_pos(1) 1.1 t3_pos(3)]);
    
    sub_4 = subplot(13,9,[78 117]);
    plot(data.DatesNum,(cumprod(1 + data.RF) - 1),'Color','b');
    set(sub_4,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Cumulative Rates');
    
    figure_title(sprintf('Market Data\nAverage Excess Return: %.2f%% | Annualized Excess Return: %.2f%%\nMAR: %.2f%% | Risk Aversion: %.2f',(data.MEAvg * 100),(data.MEAnn * 100),(data.MAR * 100),data.RiskAve));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_style_factors(date,sf)

    sf_len = width(sf);

    cols = floor(sqrt(sf_len));
    rows = floor(ceil(sf_len / cols));

    fig = figure();
    set(fig,'Name','Style Factors','Units','normalized','Position',[100 100 0.85 0.85]);

    for i = 1:sf_len
        sub = subplot(rows,cols,i);
        plot(date,sf{:,i},'Color','b');
        set(sub,'XLim',[date(1) date(end)]);
        datetick('x','yyyy','KeepLimits');
        xlabel('Time');
        ylabel(sf.Properties.VariableNames{i});
    end

    figure_title('Style Factors');

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end