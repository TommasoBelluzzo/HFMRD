% [INPUT]
% data     = The results returned by the "execute_tests" function.
% detailed = A boolean indicating whether to create a detailed plot for each time series of fund returns in the dataset (optional, default=false).

function plot_data(varargin)

    persistent ip;

    if (isempty(ip))
        ip = inputParser();
        ip.addRequired('data',@(x)validateattributes(x,{'struct'},{'nonempty'}));
        ip.addOptional('detailed',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    ip.parse(varargin{:});
    ipr = ip.Results;

    plot_data_internal(ipr.data,ipr.detailed);

end

function plot_data_internal(data,detailed)

    plot_market_data(data);
    plot_style_factors(data.DatesNum,data.StyleFactors,data.StyleFactorsNames);
    plot_funds_overview(data);
    
    if (detailed)
        benchmark = data.Benchmark;
        benchmark_cp = cumprod(1 + benchmark) - 1;

        if (numel(unique(data.Groups)) == 1)
            for i = 1:data.N
                plot_fund_details(data.DatesNum,benchmark,benchmark_cp,data.FirmReturns(:,i),data.FirmNames{i},NaN);
            end
        else
            for i = 1:data.N
                plot_fund_details(data.DatesNum,benchmark,benchmark_cp,data.FirmReturns(:,i),data.FirmNames{i},data.Groups(i));
            end
        end
    end

end

function plot_fund_details(date,benchmark,benchmark_cp,returns,firm_name,group)

    r_annualized = ((prod(1 + returns) ^ (12 / numel(returns))) - 1) * 100;
    r_avg = mean(returns) * 100;
    r_composite = cumprod(1 + returns);

    r_cumulative = r_composite - 1;
    r_cumulative_max = cummax([1; r_composite]);

    r_drawdowns = (r_composite ./ r_cumulative_max(2:end)) - 1;
    r_drawdowns_max = abs(min(r_drawdowns)) * 100;

    f = figure('Name',[firm_name ' Details'],'Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 36]);
    plot(date,returns,'Color','b');
    hold on;
        plot(date,benchmark,'Color','r');
        area_1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        area_2 = area(1,NaN,'FaceColor','b','ShowBaseLine','off');
	hold off;
    set(sub_1,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t1 = title(sub_1,'Returns','Units','normalized');
    t1_position = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_position(2) t1_position(3)]);
    
    sub_2 = subplot(13,9,[46 72]);
    plot(date,r_cumulative,'Color','b');
    hold on;
        plot(date,benchmark_cp,'Color','r');
	hold off;
    set(sub_2,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t2 = title(sub_2,'Cumulative Returns','Units','normalized');
    t2_position = get(t2,'Position');
    set(t2,'Position',[0.4783 t2_position(2) t2_position(3)]);
    
    sub_3 = subplot(13,9,[82 108]);
    plot(date,r_drawdowns,'Color','b');
    set(sub_3,'XLim',[date(1) date(end)]);
    datetick('x','yyyy','KeepLimits');
	t3 = title(sub_3,sprintf('Drawdowns (Maximum: %.2f%%)',r_drawdowns_max),'Units','normalized');
    t3_position = get(t3,'Position');
    set(t3,'Position',[0.4783 t3_position(2) t3_position(3)]);
    
    l = legend([area_1 area_2],'Benchmark',firm_name,'Location','best','Orientation','horizontal','Units','normalized');
    l_position = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_position(3)) / 2) 0.067 l_position(3) l_position(4)]);

    if (isnan(group))
        figure_title(sprintf('%s\nAverage Return: %.2f%% | Annualized Return: %.2f%%',firm_name,r_avg,r_annualized));
    else
        figure_title(sprintf('%s (Group %d)\nAverage Return: %.2f%% | Annualized Return: %.2f%%',firm_name,group,r_avg,r_annualized));
    end

    pause(0.01);
    frame = get(f,'JavaFrame');
    set(frame,'Maximized',true);

end

function plot_funds_overview(data)

    colors = [
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

    groups = data.Groups;
    groups_unique = sort(unique(groups));
    groups_len = numel(groups_unique);
    groups_legend = gobjects(groups_len,1);
    
    benchmark = data.Benchmark;
    benchmark_indices_down = benchmark <= 0;
    benchmark_down = sum(benchmark(benchmark_indices_down));
    benchmark_indices_up = benchmark > 0;
    benchmark_up = sum(benchmark(benchmark_indices_up));

    r_excess = data.FirmReturns - repmat(data.RiskFree,1,data.N);
    r_excess_annualized = (prod(1 + r_excess) .^ (12 / data.T)) - 1;
    
    capture_down = [1; zeros(data.N,1)];
    capture_up = [1; zeros(data.N,1)];
    
    for i = 1:data.N
        r = data.FirmReturns(:,i);

        capture_down(i+1) = sum(r(benchmark_indices_down)) / benchmark_down;
        capture_up(i+1) = sum(r(benchmark_indices_up)) / benchmark_up;
    end

    capture_down_x_max = (max(capture_down) + 0.3);
    capture_down_x_min = (min(capture_down) - 0.3);
    capture_up_y_max = (max(capture_up) + 0.3);
    capture_up_y_min = (min(capture_up) - 0.3);
    capture_ratio_max = max(capture_down_x_max,capture_up_y_max);
    capture_ratio_min = min(capture_down_x_min,capture_up_y_min);

    sortino_ratios = r_excess_annualized ./ (sqrt(12) .* std(data.FirmReturns));
    sortino_ratios_avg = mean(sortino_ratios);

    f = figure('Name','Funds Overview','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    b = boxplot([data.Benchmark data.FirmReturns],'Symbol','k.');
    set(findobj(sub_1,'Type','Line','Tag','Box'),'Color','k');
    set(findobj(sub_1,'Type','Line','Tag','Median'),'Color','k');
    set(findobj(sub_1,'-regexp','Tag','\w*Whisker'),'LineStyle','-');
    hold on;
        h = get(b(5,:),{'XData','YData'});
        patch(h{1,1},h{1,2},colors(1,:));
        for i = 2:size(h,1)
            patch(h{i,1},h{i,2},colors(groups(i-1)+1,:))
        end     
        sub_1_children = get(sub_1,'Children');
        set(sub_1,'Children',sub_1_children([end 1:end-1]));
        for i = 1:groups_len
            groups_legend(i) = area(0,NaN,'FaceColor',colors(i+1,:),'ShowBaseLine','off');
        end
    hold off;
    set(sub_1,'XTickLabel',[{'Benchmark'} data.FirmNames]);
    t1 = title(sub_1,'Boxes','Units','normalized');
    t1_position = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_position(2) t1_position(3)]);

    sub_2 = subplot(13,9,[73 112]);
    scatter(capture_down,capture_up,36,[colors(1,:); colors(groups+1,:)],'filled');
    hold on;
        line([capture_ratio_min capture_ratio_max],[1 1],'Color',colors(1,:));
        line([1 1],[capture_ratio_min capture_ratio_max],'Color',colors(1,:));
        text(1.03,1.07,'Benchmark','Color',colors(1,:),'Clipping','on');    
        for i = 1:data.N
            text(capture_down(i+1)+0.03,capture_up(i+1)+0.03,data.FirmNames{i},'Clipping','on')
        end
    hold off;
    grid on;
    set(sub_2,'Box','on');
    set(sub_2,'XLim',[capture_ratio_min capture_ratio_max]);
    set(sub_2,'YLim',[capture_ratio_min capture_ratio_max]);
    xlabel('Downside');
    ylabel('Upside');
    title(sub_2,'Capture Ratios');
    
    sub_3 = subplot(13,9,[78 117]);
    hold on;
        for i = 1:groups_len
            group = groups_unique(i);
            sortino_ratios_current = sortino_ratios;
            sortino_ratios_current(groups ~= group) = NaN;
            bar(1:data.N,sortino_ratios_current,'FaceColor',colors(group+1,:));
        end
        text(1:data.N,sortino_ratios,sprintfc('%.2f',sortino_ratios),'HorizontalAlignment','center','VerticalAlignment','bottom');
        line([0.5 (data.N + 0.5)],[sortino_ratios_avg sortino_ratios_avg],'Color','k');
    hold off;
    set(sub_3,'XLim',[0.5 (data.N + 0.5)],'XTick',1:data.N,'XTickLabel',data.FirmNames);
    title(sub_3,'Sharpe Ratios');

    if (groups_len > 1)
        l = legend(groups_legend,sprintfc('Group %d',groups_unique),'Location','southoutside','Orientation','horizontal');
        l_position = get(l,'Position');
        set(l,'Box','off','Position',[((1 - l_position(3)) / 2) 0.47 l_position(3) l_position(4)]);
    end

    figure_title(sprintf('Returns Overview\nObservations: %d | First: %s | Last: %s',data.T,datestr(data.DatesNum(1),'mm/yyyy'),datestr(data.DatesNum(end),'mm/yyyy')));

    pause(0.01);
    frame = get(f,'JavaFrame');
    set(frame,'Maximized',true);

end

function plot_market_data(data)

    f = figure('Name','Market Data','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[19 58]);
    plot(data.DatesNum,data.Benchmark,'Color','r');
    set(sub_1,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Returns');
	t1 = title(sub_1,sprintf('Benchmark\nAverage: %.2f%% | Annualized: %.2f%%',(data.BenchmarkAverage * 100),(data.BenchmarkAnnualized * 100)),'Units','normalized');
    t1_position = get(t1,'Position');
    set(t1,'Position',[t1_position(1) 1.1 t1_position(3)]);
    
    sub_2 = subplot(13,9,[73 112]);
    plot(data.DatesNum,(cumprod(1 + data.Benchmark) - 1),'Color','r');
    set(sub_2,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Cumulative Returns');

    sub_3 = subplot(13,9,[24 63]);
    plot(data.DatesNum,data.RiskFree,'Color','b');
    set(sub_3,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Rates');
	t3 = title(sub_3,sprintf('Risk-Free Rate\nAverage: %.2f%% | Annualized: %.2f%%',(data.RiskFreeAverage * 100),(data.RiskFreeAnnualized * 100)),'Units','normalized');
    t3_position = get(t3,'Position');
    set(t3,'Position',[t3_position(1) 1.1 t3_position(3)]);
    
    sub_4 = subplot(13,9,[78 117]);
    plot(data.DatesNum,(cumprod(1 + data.RiskFree) - 1),'Color','b');
    set(sub_4,'XLim',[data.DatesNum(1) data.DatesNum(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Cumulative Rates');
    
    figure_title(sprintf('Market Data\nAverage Excess Return: %.2f%% | Annualized Excess Return: %.2f%%\nMAR Ratio: %.2f%% | Risk Aversion: %.2f',(data.MarketExcessAverage * 100),(data.MarketExcessAnnualized * 100),(data.MAR * 100),data.RiskAversion));

    pause(0.01);
    frame = get(f,'JavaFrame');
    set(frame,'Maximized',true);

end

function plot_style_factors(date,style_factors,style_factors_names)

    style_factors_len = size(style_factors,2);

    columns = floor(sqrt(style_factors_len));
    rows = floor(ceil(style_factors_len / columns));

    f = figure('Name','Style Factors','Units','normalized','Position',[100 100 0.85 0.85]);

    for i = 1:style_factors_len
        sub = subplot(rows,columns,i);
        plot(date,style_factors(:,i),'Color','b');
        set(sub,'XLim',[date(1) date(end)]);
        datetick('x','yyyy','KeepLimits');
        title(sub,style_factors_names{i});
    end

    figure_title('Style Factors');

    pause(0.01);
    frame = get(f,'JavaFrame');
    set(frame,'Maximized',true);

end