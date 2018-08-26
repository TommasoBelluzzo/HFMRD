% [INPUT]
% td = A cell array of structures returned by the "execute_tests" function when the "tab" parameter is set to "false".

function plot_results(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('td',@(x)validateattributes(x,{'cell'},{'nonempty'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    td = res.td;

    plot_results_internal(td);

end

function t = figure_title(str)

    fig = gcf();
    fig_fts = get(fig,'DefaultAxesFontSize') + 4;
    fig_uni = get(fig,'Units');
    
    if (~strcmp(fig_uni,'pixels'))
        set(fig,'Units','pixels');
        fig_pos = get(fig,'Position');
        set(fig,'Units',fig_uni);
    else
        fig_pos = get(fig,'Position');
    end

    ff = ((fig_fts - 4) * 6.35) / fig_pos(4);

    tit = NaN;
    y_max = 0;
    y_min = 1;

    h = findobj(fig,'Type','axes');
    h_len = length(h);
    h_pos = zeros(h_len,4);

    for i = 1:h_len
        h_cur = h(i);
        
        fig_pos = get(h_cur,'Position');
        h_pos(i,:) = fig_pos;

        if (~strcmp(get(h_cur,'Tag'),'suptitle'))
            fig_y = fig_pos(2);
            fig_hei = fig_pos(4);
            
            if (fig_y < y_min)
                y_min = fig_y - (ff / 15);
            end

            if ((fig_hei + fig_y) > y_max)
                y_max = fig_hei + fig_y + (ff / 10);
            end
        else
            tit = h_cur;
        end
    end

    if (y_max > 0.92)
        scl = (0.92 - y_min) / (y_max - y_min);

        for i = 1:h_len
            fig_pos = h_pos(i,:);
            fig_pos(2) = ((fig_pos(2) - y_min) * scl) + y_min;
            fig_pos(4) = (fig_pos(4) * scl) - ((1 - scl) * (ff / 15));

            set(h(i),'Position',fig_pos);
        end
    end

    np = get(fig,'NextPlot');
    set(fig,'NextPlot','add');

    if (ishghandle(tit))
        delete(tit);
    end

    axes('Position',[0 1 1 1],'Tag','suptitle','Visible','off');
    t_int = text(0.50,-0.05,str,'HorizontalAlignment','center','FontSize',fig_fts);

    set(fig,'NextPlot',np);

    axes(gca());

    if (nargout == 1)
        t = t_int;
    end

end

function plot_results_internal(td)

    fig = findobj(allchild(0),'flat','Tag','TestResults');
    
    if (~isempty(fig))
        close(fig);
    end

    n = size(td,2);
    td1 = td{1,1};
    
    heads = cell(n,1);
    fail = ones(7,n) .* repmat((7:-1:1).',1,3);
    succ = ones(7,n) .* repmat((7:-1:1).',1,3);
    
    for i = 1:n
        coef = 0;
        
        for j = 1:7
            t = td{j,i};
            
            coef = coef + t.FailCoef;
            
            if (t.Fail)
                succ(j,i) = NaN;
            else
                fail(j,i) = NaN;
            end
        end

        heads{i} = sprintf('%s\nAnomaly Score: %.1f/7 (%.0f%%)',t.Frm,coef,round(((coef / 7) * 100)));
    end

    fig = figure();
    set(fig,'Name','Test Results','Tag','TestResults','Units','normalized','Position',[100 100 0.85 0.85]);

    hold on;
        for i = 1:n
            for j = 1:7
                b1 = bar(i,fail(j,i),'FaceColor',[0.800 0.000 0.000]);
                b2 = bar(i,succ(j,i),'FaceColor',[0.298 0.800 0.000]);
                set([b1 b2],'BarWidth',1,'ButtonDownFcn',{@plot_results_switch,td{j,i}});
            end
        end
    hold off;

    a = gca();
    set(a,'FontSize',12);
    set(a,'XLim',[0.5 (n + 0.5)],'XTick',1:n,'XTickLabel',[]);
    set(a,'YLim',[0 7],'YTick',0.5:1:7,'YTickLabel',{'Data Quality' 'Digits Conformity' 'Discontinuity At Zero' 'December Spike' 'Bias Ratio' 'Serial Correlation' 'Low Correlation'});

    for i = 1:n
        text(i,-0.2, heads{i},'FontSize',12,'HorizontalAlignment','center');
    end
    
    figure_title(sprintf('Test Results\na: %.2f | Observations: %d | Simulations: %d',td1.Par.A,td1.Par.Obs,td1.Par.Sims));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_results_switch(obj,evd,td) %#ok<INUSL>

    fig = findobj(allchild(0),'flat','Tag','TestPlot');
    
    if (~isempty(fig))
        close(fig);
    end

    switch (td.Type)          
        case 'BR'
            plot_test_bias_ratio(td);
        case 'DC'
            plot_test_digits_conformity(td);
        case 'DQ'
            plot_test_data_quality(td);
        case 'DS'
            plot_test_december_spike(td);
        case 'DZ'
            plot_test_discontinuity_at_zero(td);
        case 'LC'
            plot_test_low_conformity(td);
        case 'SC'
            plot_test_serial_correlation(td);
        otherwise
            error('Unrecognized test type ''%s''.',test);
    end

end

function plot_test_bias_ratio(td)

    data = td.Data;
    
    [bins,edgs] = histcounts(data.H0,10);
    edgs_wid = edgs(2) - edgs(1);
    edgs_seq = [fliplr((edgs(1)-edgs_wid):-edgs_wid:0) edgs(1):edgs_wid:max(edgs(end),((ceil(data.BR * 2) / 2) * 1.05))];

    all_tit = sprintf('Bias Ratio Test on %s',td.Frm);

    if (td.Fail)
        all_clr = [0.800 0.000 0.000];
        all_res = 'Failure';
    else
        all_clr = [0.298 0.800 0.000];
        all_res = 'Success';
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 117]);
    histogram(sub_1,'BinEdges',edgs,'BinCounts',bins,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        ymax = get(sub_1,'YLim');
        line([data.BR data.BR],[0 ymax(2)],'Color',all_clr,'LineWidth',1.5);
    hold off;
    set(sub_1,'XLim',[edgs_seq(1) edgs_seq(end)],'XTick',edgs_seq,'XTickLabel',sprintfc('%.2f',edgs_seq));
    xlabel('Simulated Bias Ratios');
    ylabel('Counts');
    
    st = figure_title(sprintf('%s\nValue: %.2f | p-Value: %.4f\nResult: %s',all_tit,data.BR,data.PVal,all_res));
    st_pos = get(st,'Position');
    set(st,'Color',all_clr,'Position',[0.5170 -0.080 st_pos(3)]);

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_data_quality(td)

    data = td.Data;
    
    all_coef = td.Flags * td.FailCoef;
    all_flag = td.Flags;
    all_tit = sprintf('Data Quality Test on %s',td.Frm);
    
    if (all_coef == 1)
        all_plur = '';
    else
        all_plur = 's';
    end
    
    if (td.Fail)
        all_res = 'Failure';
    else
        all_res = 'Success';
    end

    neg_len = numel(data.NegProbX);
    neg_seq = 1:neg_len;
    neg_int = linspace(1,neg_len,1000);
    neg_int_x = interp1(neg_seq,data.NegProbX,neg_int,'spline');
    neg_int_f = interp1(neg_seq,data.NegProbF,neg_int,'spline');
    neg_curr_x = data.NegProbX(1:data.NegVI);
    neg_curr_f = data.NegProbF(1:data.NegVI);
    neg_len = numel(neg_curr_x);
    neg_seq = 1:neg_len;
    neg_int = linspace(1,neg_len,1000);
    neg_curr_int_x = interp1(neg_seq,neg_curr_x,neg_int,'spline');
    neg_curr_int_f = interp1(neg_seq,neg_curr_f,neg_int,'spline');
	neg_txt = sprintf('Negative Returns\nCount: %d | Critical Value: %d',data.NegVal,data.NegTest);

    if (data.NegFail)
        neg_clr = [0.800 0.000 0.000];
    else
        neg_clr = [0.298 0.800 0.000];
    end

    zer_len = numel(data.ZerProbX);
    zer_seq = 1:zer_len;
    zer_int = linspace(1,zer_len,1000);
    zer_int_x = interp1(zer_seq,data.ZerProbX,zer_int,'spline');
    zer_int_f = interp1(zer_seq,data.ZerProbF,zer_int,'spline');
    zer_curr_x = data.ZerProbX(1:data.ZerVI);
    zer_curr_f = data.ZerProbF(1:data.ZerVI);
    zer_len = numel(zer_curr_x);
    zer_seq = 1:zer_len;
    zer_int = linspace(1,zer_len,1000);
    zer_curr_int_x = interp1(zer_seq,zer_curr_x,zer_int,'spline');
    zer_curr_int_f = interp1(zer_seq,zer_curr_f,zer_int,'spline');
    zer_txt = sprintf('Zero Returns\nCount: %d | Critical Value: %d',data.ZerVal,data.ZerTest);
    
    if (data.ZerFail)
        zer_clr = [0.800 0.000 0.000];
    else
        zer_clr = [0.298 0.800 0.000];
    end
    
    pai_txt = sprintf('Return Pairs\nCount: %d | Critical Value: %d',data.PaiVal,data.PaiTest);
    
    if (data.PaiFail)
        pai_clr = [0.800 0.000 0.000];
    else
        pai_clr = [0.298 0.800 0.000];
    end
    
    str_txt = sprintf('Return Strings\nCount: %d | Critical Value: %d',data.StrVal,data.StrTest);
    
    if (data.StrFail)
        str_clr = [0.800 0.000 0.000];
    else
        str_clr = [0.298 0.800 0.000];
    end
    
    uni_txt = sprintf('Unique Returns\nCount: %d | Critical Value: %d',data.UniVal,data.UniTest);
    
    if (data.UniFail)
        uni_clr = [0.800 0.000 0.000];
    else
        uni_clr = [0.298 0.800 0.000];
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 49]);
    plot(neg_int_x,neg_int_f,'Color',[0.239 0.149 0.659]);
    hold on;
        x_lim = get(sub_1,'XLim');
        plot(neg_curr_int_x,neg_curr_int_f,'Color',neg_clr);
        line([neg_int_x(1) x_lim(2)],[data.NegCV data.NegCV],'Color',[0.800 0.450 0.000]);
        l1 = area(0,NaN,'FaceColor',[0.800 0.450 0.000],'ShowBaseLine','off');
        l2 = area(0,NaN,'FaceColor',[0.800 0.000 0.000],'ShowBaseLine','off');
        l3 = area(0,NaN,'FaceColor',[0.298 0.800 0.000],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[neg_int_x(1) neg_int_x(end)]);
    set(sub_1,'YLim',[-0.025 1.025]);
    grid on;
    xlabel('Number of Negative Returns');
    ylabel('Simulated Cumulative Probability');
    title(sub_1,neg_txt,'Color',neg_clr);
    
    sub_2 = subplot(13,9,[15 54]);
    plot(zer_int_x,zer_int_f,'Color',[0.239 0.149 0.659]);
    hold on;
        x_lim = get(sub_2,'XLim');
        plot(zer_curr_int_x,zer_curr_int_f,'Color',zer_clr);
        line([zer_int_x(1) x_lim(2)],[data.ZerCV data.ZerCV],'Color',[0.800 0.450 0.000]);
    hold off;
    set(sub_2,'XGrid','on');
    set(sub_2,'XLim',[zer_int_x(1) zer_int_x(end)]);
    set(sub_2,'YLim',[-0.025 1.025]);
    grid on;
    xlabel('Number of Zero Returns');
    ylabel('Simulated Cumulative Probability');
    title(sub_2,zer_txt,'Color',zer_clr);

    sub_3 = subplot(13,9,[73 101.3]);
    histogram(sub_3,data.PaiH0,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        y_lim = get(sub_3,'YLim');
        line([data.PaiTest data.PaiTest],[0 y_lim(2)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
        line([data.PaiVal data.PaiVal],[0 y_lim(2)],'Color',pai_clr,'LineWidth',1.5);
    hold off;
    xlabel('Simulated Return Pairs');
    ylabel('Counts');
    title(sub_3,pai_txt,'Color',pai_clr);
    
    sub_4 = subplot(13,9,[75.7 104]);
    histogram(sub_4,data.StrH0,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        y_lim = get(sub_4,'YLim');
        line([data.StrTest data.StrTest],[0 y_lim(2)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
        line([data.StrVal data.StrVal],[0 y_lim(2)],'Color',str_clr,'LineWidth',1.5);
    hold off;
    xlabel('Simulated Return Strings');
    ylabel('Counts');
    title(sub_4,str_txt,'Color',str_clr);
    
    sub_5 = subplot(13,9,[79 108]);
    histogram(sub_5,data.UniH0,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        y_lim = get(sub_5,'YLim');
        line([data.UniTest data.UniTest],[0 y_lim(2)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
        line([data.UniVal data.UniVal],[0 y_lim(2)],'Color',uni_clr,'LineWidth',1.5);
    hold off;
    xlabel('Simulated Unique Returns');
    ylabel('Counts');
    title(sub_5,uni_txt,'Color',uni_clr);

    l = legend([l1 l2 l3],sprintf('Critical Values (%.1f%%)',((1 - td.Par.A) * 100)),'Failed Flags','Successful Flags','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('%s\nResult: %s (a: %.2f, %d raised flag%s out of %d)',all_tit,all_res,td.Par.A,all_coef,all_plur,all_flag));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_december_spike(td)

    data = td.Data;

    all_tit = sprintf('December Spike Test on %s',td.Frm);

    if (td.Fail)
        all_clr = [0.800 0.000 0.000];
        all_res = 'Failure';
    else
        all_clr = [0.298 0.800 0.000];
        all_res = 'Success';
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 108]);
    histogram(sub_1,data.H0Spr,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;   
        y_lim = get(sub_1,'YLim');
        line([data.H0Prc data.H0Prc],[0 y_lim(2)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
        line([data.FrmSpr data.FrmSpr],[0 y_lim(2)],'Color',all_clr,'LineWidth',1.5);
        l1 = area(0,NaN,'FaceColor',[0.800 0.450 0.000],'ShowBaseLine','off');
        l2 = area(0,NaN,'FaceColor',all_clr,'ShowBaseLine','off');
    hold off;
    xlabel('Simulated Spreads');
    ylabel('Counts');
    
    l = legend([l1 l2],sprintf('Critical Value (%.1f%%)',((1 - td.Par.A) * 100)),'Spread','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);
    
    st = figure_title(sprintf('%s\nCount(M=12): %d (%.2f%%) | Mean(M=12): %.4f | Mean(M<>12): %.4f | Spread: %.4f\nResult: %s',all_tit,data.Dec,(data.DecPrc * 100),data.FrmAvgDec,data.FrmAvgOth,data.FrmSpr,all_res));
    st_pos = get(st,'Position');
    set(st,'Color',all_clr,'Position',[0.5170 -0.080 st_pos(3)]);

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_digits_conformity(td)

    data = td.Data;
    
    all_coef = td.Flags * td.FailCoef;
    all_flag = td.Flags;
    all_tit = sprintf('Digits Conformity Test on %s',td.Frm);
    
    if (all_coef == 1)
        all_plur = '';
    else
        all_plur = 's';
    end
    
    if (td.Fail)
        all_res = 'Failure';
    else
        all_res = 'Success';
    end
    
    fst_txt = sprintf('Benford''s Law Conformity of First Digits\nX2: %.4f | p-Value: %.4f',data.FstChi2,data.FstPVal);
    
    if (data.FstFail)
        fst_clr = [0.800 0.000 0.000];
    else
        fst_clr = [0.298 0.800 0.000];
    end
    
    lst_txt = sprintf('Uniform Distribution Conformity of Last Digits\nX2: %.4f | p-Value: %.4f',data.LstChi2,data.LstPVal);
    
    if (data.LstFail)
        lst_clr = [0.800 0.000 0.000];
    else
        lst_clr = [0.298 0.800 0.000];
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 103]);
    set(bar(data.FstEmpP),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        stem(data.FstTheP,'Color',[0.800 0.450 0.000],'Marker','.','MarkerSize',20);
        l1 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.800 0.450 0.000],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[0.5 9.5],'XTick',1:9,'XTickLabel',sprintfc('%d',1:9));
    set(sub_1,'YLim',[0 (max(max(data.FstEmpP),max(data.FstTheP)) * 1.05)]);
    xlabel('Digits');
    ylabel('Frequencies');
    title(sub_1,fst_txt,'Color',fst_clr);

    sub_2 = subplot(13,9,[15 108]);
    set(bar(data.LstEmpP),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        line([0.5 10.5],[data.LstTheP(1) data.LstTheP(1)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    hold off;
    set(sub_2,'XLim',[0.5 10.5],'XTick',1:10,'XTickLabel',sprintfc('%d',0:9));
    set(sub_2,'YLim',[0 (max(max(data.LstEmpP),max(data.LstTheP)) * 1.05)]);
    xlabel('Digits');
    ylabel('Frequencies');
    title(sub_2,lst_txt,'Color',lst_clr);
    
    l = legend([l1 l2],'Empirical Frequencies','Theorical Frequencies','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('%s\nResult: %s (a: %.2f, %d raised flag%s out of %d)',all_tit,all_res,td.Par.A,all_coef,all_plur,all_flag));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_discontinuity_at_zero(td)

    data = td.Data;
    
    bins = data.Bins;
    idxz = data.IdxZ;
    
    bins_hl = zeros(1,numel(bins));
    bins_hl(idxz - 1) = bins(idxz - 1);
    bins_hl(idxz) = bins(idxz);

    all_tit = sprintf('Discontinuity At Zero Test on %s',td.Frm);

    if (td.Fail)
        all_clr = [0.800 0.000 0.000];
        all_res = 'Failure';
    else
        all_clr = [0.298 0.800 0.000];
        all_res = 'Success';
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 117]);
    histogram(sub_1,'BinEdges',data.Edges,'BinCounts',bins,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        histogram(sub_1,'BinEdges',data.Edges,'BinCounts',bins_hl,'FaceAlpha',1,'FaceColor',all_clr);
    hold off;
    set(sub_1,'XLim',[-0.03 0.03]);
    xlabel('Returns [-0.03 +0.03]');
    ylabel('Counts');
    
    st = figure_title(sprintf('%s\nDelta: %.1f | Delta Variance: %.4f | p-Value: %.4f\nResult: %s',all_tit,data.Diff,data.DiffVar,data.PVal,all_res));
    st_pos = get(st,'Position');
    set(st,'Color',all_clr,'Position',[0.5170 -0.080 st_pos(3)]);

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_low_conformity(td)

    data = td.Data;
    
    all_coef = td.Flags * td.FailCoef;
    all_flag = td.Flags;
    all_tit = sprintf('Low Correlation Test on %s',td.Frm);
    
    if (all_coef == 1)
        all_plur = '';
    else
        all_plur = 's';
    end
    
    if (td.Fail)
        all_res = 'Failure';
    else
        all_res = 'Success';
    end
    
    idx_mdl = data.IdxMod;
    idx_a = idx_mdl.Coefficients{1,1};
    idx_b = idx_mdl.Coefficients{2,1};
    idx_pval = idx_mdl.Coefficients{2,4};
    idx_r2a = idx_mdl.Rsquared.Adjusted * 100;
    idx_txt = sprintf('Equally Weighted Index Flag\nModel: Rt = %.4f + %.4f * EWIt\nAdjusted R2: %.2f%% | B:  %.4f | B p-Value: %.4f',idx_a,idx_b,idx_r2a,idx_b,idx_pval);
    
    if (data.IdxFail)
        idx_clr = [0.800 0.000 0.000];
    else
        idx_clr = [0.298 0.800 0.000];
    end
    
    max_mdl = data.MaxMod;

    if (isempty(data.MaxCP))
        if (numel(data.MaxComb) == 1)
            max_sf_txt = ['Maximum Style Factor: ' data.MaxComb{1}];
        else
            max_sf_txt = sprintf('Maximum Style Factors: %s',strjoin(data.MaxComb,', '));
        end
    else
        max_coms_len = (numel(data.MaxComb) - 1) / 2;
        
        max_sf_txt = ['Change Point: ' datestr(data.MaxCP,'mm/yyyy')];
        
        if (max_coms_len == 1)
            max_sf_txt = [max_sf_txt ' | Maximum Style Factors: ' strjoin(data.MaxComb,' ')];
        else
            max_sf_txt = [max_sf_txt ' | Maximum Style Factors: ' strjoin(data.MaxComb(1:max_coms_len),', ') ' ' data.MaxComb{max_coms_len+1} ' ' strjoin(data.MaxComb((max_coms_len + 2):end),', ')];
        end
    end

    max_txt = sprintf('Maximum R2 Flag\n%s\nMaximum Adjusted R2: %.2f%% | Critical Value: %.2f%%',max_sf_txt,(data.MaxVal * 100),(data.MaxTest * 100));
    
    if (data.MaxFail)
        max_clr = [0.800 0.000 0.000];
    else
        max_clr = [0.298 0.800 0.000];
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[19 49]);
    p1 = plot(idx_mdl);
    set(p1(1),'Color',[0.239 0.149 0.659],'MarkerSize',10);
    set(p1(2),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
	set(p1(3:end),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    delete(findobj(fig,'Type','Legend'));
    xlabel('X');
    ylabel('Y');
	t1 = title(sub_1,idx_txt,'Color',idx_clr,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[t1_pos(1) 1.1 t1_pos(3)]);
    legend([p1(1) p1(2) p1(3)],'Data','Fit','Confidence Bounds (95%)','Location','best','Orientation','horizontal','Units','normalized');

    sub_2 = subplot(13,9,[64 103]);
    p_idx = plot(td.Par.Dates,data.IdxMod.Variables.y,'Color',[0.239 0.149 0.659]);
    hold on;
        p_ewi = plot(td.Par.Dates,data.IdxEWI,'Color',[0.800 0.450 0.000]);
    hold off;
    set(sub_2,'XLim',[td.Par.Dates(1) td.Par.Dates(end)]);
    datetick('x','yyyy','KeepLimits');
    xlabel('Time');
    ylabel('Returns');
    legend([p_idx p_ewi],td.Frm,'Equally Weighted Index','Location','best','Orientation','vertical','Units','normalized');
    
    sub_3 = subplot(13,9,[24 54]);
    p2 = plot(max_mdl);
    set(p2(1),'Color',[0.239 0.149 0.659],'MarkerSize',10);
    set(p2(2),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    set(p2(3:end),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    xlabel('X');
    ylabel('Y');
    t3 = title(sub_3,max_txt,'Color',max_clr,'Units','normalized');
    t3_pos = get(t3,'Position');
    set(t3,'Position',[t3_pos(1) 1.1 t3_pos(3)]);
    legend([p2(1) p2(2) p2(3)],'Data','Fit','Confidence Bounds (95%)','Location','best','Orientation','horizontal','Units','normalized');

    sub_4 = subplot(13,9,[69 108]);
    histogram(sub_4,data.MaxH0,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659]);
    hold on;
        y_lim = get(sub_4,'YLim');
        line([data.MaxTest data.MaxTest],[0 y_lim(2)],'Color',[0.800 0.450 0.000],'LineWidth',1.5);
        line([data.MaxVal data.MaxVal],[0 y_lim(2)],'Color',max_clr,'LineWidth',1.5);
    hold off;
    xlabel('Simulated Maximum Adjusted R2');
    ylabel('Counts');

    figure_title(sprintf('%s\nResult: %s (a: %.2f, %d raised flag%s out of %d)',all_tit,all_res,td.Par.A,all_coef,all_plur,all_flag));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_test_serial_correlation(td)

    data = td.Data;
    
    all_coef = td.Flags * td.FailCoef;
    all_flag = td.Flags;
    all_tit = sprintf('Serial Correlation Test on %s',td.Frm);
    
    if (all_coef == 1)
        all_plur = '';
    else
        all_plur = 's';
    end
    
    if (td.Fail)
        all_res = 'Failure';
    else
        all_res = 'Success';
    end
    
    unc_mdl = data.UncMod;
    unc_a = unc_mdl.Coefficients{1,1};
    unc_b = unc_mdl.Coefficients{2,1};
    unc_pval = unc_mdl.Coefficients{2,4};
    unc_r2a = unc_mdl.Rsquared.Adjusted * 100;
    
    if (sign(unc_b) >= 0)
        sign_b = '+';
    else
        sign_b = '-';
    end
    
    unc_txt = sprintf('Unconditional Serial Correlation\nModel: Rt = %.4f %s %.4f * R(t-1)\nAdjusted R2: %.2f%% | B: %.4f | B p-Value: %.4f',unc_a,sign_b,abs(unc_b),unc_r2a,unc_b,unc_pval);
    
    if (data.UncFail)
        unc_clr = [0.800 0.000 0.000];
    else
        unc_clr = [0.298 0.800 0.000];
    end
    
    con_mdl = data.ConMod;
    con_a = con_mdl.Coefficients{1,1};
    con_b_pos = con_mdl.Coefficients{2,1};
    con_b_neg = con_mdl.Coefficients{3,1};
    con_pval = con_mdl.Coefficients{3,4};
    con_r2a = con_mdl.Rsquared.Adjusted * 100;
    
    if (sign(con_b_pos) >= 0)
        sign_b_pos = '+';
    else
        sign_b_pos = '-';
    end
    
    if (sign(con_b_neg) >= 0)
        sign_b_neg = '+';
    else
        sign_b_neg = '-';
    end
    
    con_txt = sprintf('Conditional Serial Correlation\nModel: Rt = %.4f %s %.4f * R(t-1) %s %.4f * [1 - I(t-1)] * R(t-1)\nAdjusted R2: %.2f%% | B2: %.4f | B2 p-Value: %.4f',con_a,sign_b_pos,abs(con_b_pos),sign_b_neg,abs(con_b_neg),con_r2a,con_b_neg,con_pval);
    
    if (data.ConFail)
        con_clr = [0.800 0.000 0.000];
    else
        con_clr = [0.298 0.800 0.000];
    end

    fig = figure();
    set(fig,'Name',all_tit,'Tag','TestPlot','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[19 103]);
    p1 = plot(unc_mdl);
    set(p1(1),'Color',[0.239 0.149 0.659],'MarkerSize',10);
    set(p1(2),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
	set(p1(3:end),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    delete(findobj(fig,'Type','Legend'));
    xlabel('X');
    ylabel('Y');
	title(sub_1,unc_txt,'Color',unc_clr);

    sub_2 = subplot(13,9,[24 108]);
    p2 = plot(con_mdl);
    set(p2(1),'Color',[0.239 0.149 0.659],'MarkerSize',10);
    set(p2(2),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
	set(p2(3:end),'Color',[0.800 0.450 0.000],'LineWidth',1.5);
    delete(findobj(fig,'Type','Legend'));
    xlabel('X');
    ylabel('Y');
	title(sub_2,con_txt,'Color',con_clr);
    
    l = legend([p1(1) p1(2) p1(3)],'Data','Fit','Confidence Bounds (95%)','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('%s\nResult: %s (a: %.2f, %d raised flag%s out of %d)',all_tit,all_res,td.Par.A,all_coef,all_plur,all_flag));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end