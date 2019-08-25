% [INPUT]
% data              = A structure representing the dataset.
% a                 = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.10).
% simulations       = An integer (>= 1000) representing the number of Monte Carlo simulations to perform (optional, default=10000).
% style_factors_max = An integer [1,3] representing the maximum number of style factors that can be used simultaneously in regression models (optional, default=3).
% r2_switch         = A boolean indicating whether to switch style factors during the calculation of maximum R2 in the low correlation test (optional, default=false).
% decimals          = An integer [2,6] representing the number of decimal places to consider when performing the digit tests (optional, default=4).
%
% [OUTPUT]
% results           = A cell array of structures with the quantitative details of each test.
% summary           = A table that displays the outcome and the anomaly score of each test.

function [results,summary] = execute_tests(varargin)

    persistent ip;

    if (isempty(ip))
        ip = inputParser();
        ip.addRequired('data',@(x)validateattributes(x,{'struct'},{'nonempty'}));
        ip.addOptional('a',0.10,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        ip.addOptional('simulations',10000,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1000}));
        ip.addOptional('style_factors_max',3,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1,'<=',6}));
        ip.addOptional('r2_switch',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
        ip.addOptional('decimals',4,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',2,'<=',6}));
    end
    
    ip.parse(varargin{:});

    ipr = ip.Results;
    data = ipr.data;
    style_factors_max = max([ipr.style_factors_max size(data.StyleFactors,2)]);
    
    nargoutchk(1,2);
    
    if (nargout == 1)
        results = execute_tests_internal(data,ipr.a,ipr.simulations,style_factors_max,ipr.r2_switch,ipr.decimals);
    else
        [results,summary] = execute_tests_internal(data,ipr.a,ipr.simulations,style_factors_max,ipr.r2_switch,ipr.decimals);
    end

end

function [results,summary] = execute_tests_internal(data,a,simulations,style_factors_max,r2_switch,decimals)

    params = struct();
    params.A = a;
    params.Dates = data.DatesNum;
    params.Obs = data.Obs;
    params.Sims = simulations;

    style_factors_comb = nchoosek(1:size(data.StyleFactors,2),style_factors_max);

    output = cell(7,data.Frms);
    
    for i = 1:data.Frms
        firm_name = data.FirmNames{i};
        
        returns = data.FrmsRet(:,i);
        returns_h0 = round(normrnd(mean(returns),std(returns),[data.Obs simulations]),4);
        returns_other = data.FrmsRet(:,(~strcmp(data.FirmNames,data.FirmNames(i)) & (data.Grps == data.Grps(i))));

        [output_lc,returns_fitted] = test_low_correlation(params,firm_name,returns,returns_h0,returns_other,data.StyleFactors,style_factors_comb,data.StyleFactorsNames,r2_switch);

        output{1,i} = output_lc;
        output{2,i} = test_serial_correlation(params,firm_name,returns,returns_fitted);
        output{3,i} = test_bias_ratio(params,firm_name,returns,returns_h0);
        output{4,i} = test_december_spike(params,firm_name,returns,returns_h0);
        output{5,i} = test_discontinuity_at_zero(params,firm_name,returns);
        output{6,i} = test_digits_conformity(params,firm_name,returns,decimals);
        output{7,i} = test_data_quality(params,firm_name,returns,returns_h0);
    end

    results = output;

    if (nargout == 2)
        summary = cell(data.Frms,1);

        for j = 1:data.Frms
            coefficients = zeros(8,1);
            failures = false(8,1);

            for i = 1:7
                test_curr = output{i,j};
                coefficients(i) = test_curr.FailCoef;
                failures(i) = test_curr.Fail;
            end

            coefficients(8) = sum(coefficients(1:7));
            failures(8) = coefficients(8) > 3.5;

            summary{j} = table(failures,coefficients,'VariableNames',{'Failure' 'Coefficient'});
        end

        summary = table(summary{:},'VariableNames',data.FirmNames);
        summary.Properties.RowNames = {'Low Correlation' 'Serial Correlation' 'Bias Ratio' 'December Spike' 'Discontinuity At Zero' 'Digits Conformity' 'Data Quality' 'Total'};
    end

end

function output = test_bias_ratio(params,firm_name,returns,returns_h0)

    upper_bound = std(returns);
    lower_bound = -upper_bound;
    bias_ratio = sum((returns >= 0) & (returns <= upper_bound)) / sum((returns >= lower_bound) & (returns < 0));

    bias_ratios_h0 = zeros(size(returns_h0,2),1);
    
    for i = 1:size(returns_h0,2)
        r = returns_h0(:,i);
        
        upper_bound = std(r);
        lower_bound = -upper_bound;
        bias_ratios_h0(i) = sum((r >= 0) & (r <= upper_bound)) / sum((r >= lower_bound) & (r < 0));
    end
    
    pval = sum(bias_ratios_h0 >= (bias_ratio - 1e-8)) / size(returns_h0,2);
    pval = max([0 min([pval 1])]);
    
    output = struct();

    output.Type = 'Bias Ratio';
    output.Flags = 1;
    output.Par = params;
    output.Frm = firm_name;

    output.Fail = (bias_ratio >= 2.5) & (pval < params.A);
    output.FailCoef = output.Fail;
    
    output.Data = struct();
    output.Data.BR = bias_ratio;
    output.Data.H0 = bias_ratios_h0;
    output.Data.PVal = pval;

end

function res = test_data_quality(par,frm,ret,ret_h0)

    a = par.A;
    n = par.Obs;

    ret = round(ret,2);
    ret_diff = diff([0; find(diff(ret)); n]);

    ret_h0 = round(ret_h0,2);
    ret_h0_len = size(ret_h0,2);

    ret_neg = sum(ret < 0);
    ret_h0_neg = sum(ret_h0 < 0);
    ret_h0_neg_prc = prctile(ret_h0_neg,(a * 100));

    [ret_h0_neg_f,ret_h0_neg_x] = ecdf(ret_h0_neg);
    ret_h0_neg_x(1) = max([0 (ret_h0_neg_x(1) - 1)]);
    ret_h0_neg_cv = ret_h0_neg_f(ret_h0_neg_x == ret_h0_neg_prc);
    ret_h0_neg_vi = find(ret_h0_neg_x == ret_neg);

    ret_pai = sum(ret_diff(ret_diff > 1) - 1);
    ret_h0_pai = zeros(ret_h0_len,1);
    
    for i = 1:ret_h0_len
        ret_h0_d = diff([0; find(diff(ret_h0(:,i))); n]);
        ret_h0_pai(i) = sum(ret_h0_d(ret_h0_d > 1) - 1);
    end
    
    ret_h0_pai_prc = prctile(ret_h0_pai,((1 - a) * 100));

    ret_str = max([0; ret_diff(ret_diff > 1)]);
    ret_h0_str = zeros(ret_h0_len,1);
    
    for i = 1:ret_h0_len
        ret_h0_d = diff([0; find(diff(ret_h0(:,i))); n]);
        ret_h0_str(i) = max([0; ret_h0_d(ret_h0_d > 1)]);
    end

    ret_h0_str_prc = prctile(ret_h0_str,((1 - a) * 100));

    ret_uni = numel(unique(ret));
    ret_h0_uni = zeros(ret_h0_len,1);
    
    for i = 1:ret_h0_len
        ret_h0_uni(i) = numel(unique(ret_h0(:,i)));
    end
    
    ret_h0_uni_prc = prctile(ret_h0_uni,(a * 100));

    ret_zer = sum(ret == 0);
    ret_h0_zer = sum(ret_h0 == 0);
    ret_h0_zer_prc = prctile(ret_h0_zer,((1 - a) * 100));

    [ret_h0_zer_f,ret_h0_zer_x] = ecdf(ret_h0_zer);
    ret_h0_zer_f = 1 - ret_h0_zer_f;
    ret_h0_zer_x(1) = max([0 (ret_h0_zer_x(1) - 1)]);
    ret_h0_zer_cv = ret_h0_zer_f(ret_h0_zer_x == ret_h0_zer_prc);
    ret_h0_zer_vi = find(ret_h0_zer_x == ret_zer);
    
    fail_neg = (ret_neg < ret_h0_neg_prc);
    fail_pai = (ret_pai > ret_h0_pai_prc);
    fail_str = (ret_str > ret_h0_str_prc);
    fail_uni = (ret_uni < ret_h0_uni_prc);
    fail_zer = (ret_zer > ret_h0_zer_prc);

    res = struct();
    res.Type = 'DQ';
    res.Flags = 5;
    res.Par = par;
    res.Frm = frm;
    res.Fail = fail_neg | fail_pai | fail_str | fail_uni | fail_zer;
    res.FailCoef = (fail_neg + fail_pai + fail_str + fail_uni + fail_zer) / 5;
    
    res.Data = struct();
    res.Data.NegFail = fail_neg;
    res.Data.NegCV = ret_h0_neg_cv;
    res.Data.NegH0 = ret_h0_neg;
    res.Data.NegProbF = ret_h0_neg_f;
    res.Data.NegProbX = ret_h0_neg_x;
    res.Data.NegTest = ret_h0_neg_prc;
    res.Data.NegVal = ret_neg;
    res.Data.NegVI = ret_h0_neg_vi;
    res.Data.PaiFail = fail_pai;
    res.Data.PaiH0 = ret_h0_pai;
    res.Data.PaiTest = ret_h0_pai_prc;
    res.Data.PaiVal = ret_pai;
    res.Data.StrFail = fail_str;
    res.Data.StrH0 = ret_h0_str;
    res.Data.StrTest = ret_h0_str_prc;
    res.Data.StrVal = ret_str;
    res.Data.UniFail = fail_uni;
    res.Data.UniH0 = ret_h0_uni;
    res.Data.UniTest = ret_h0_uni_prc;
    res.Data.UniVal = ret_uni;
    res.Data.ZerFail = fail_zer;
    res.Data.ZerCV = ret_h0_zer_cv;
    res.Data.ZerH0 = ret_h0_zer;
    res.Data.ZerProbF = ret_h0_zer_f;
    res.Data.ZerProbX = ret_h0_zer_x;
    res.Data.ZerTest = ret_h0_zer_prc;
    res.Data.ZerVal = ret_zer;
    res.Data.ZerVI = ret_h0_zer_vi;

end

function output = test_december_spike(params,firm_name,returns,returns_h0)

    indices = (month(params.Dates) == 12);
    
    returns_december = returns(indices,:);
    returns_december_avg = mean(returns_december);
    returns_other = returns(~indices,:);
    returns_other_avg = mean(returns_other);
    spread = returns_december_avg - returns_other_avg;
    
    returns_h0_december = returns_h0(indices,:);
    returns_h0_december_avg = mean(returns_h0_december);
    returns_h0_other = returns_h0(~indices,:);
    returns_h0_other_avg = mean(returns_h0_other);
    spread_h0 = returns_h0_december_avg - returns_h0_other_avg;
    percentile_h0 = prctile(spread_h0,((1 - params.A) * 100));
    
    output = struct();

    output.Type = 'December Spike';
    output.Flags = 1;
    output.Par = params;
    output.Frm = firm_name;

    output.Fail = (spread > percentile_h0);
    output.FailCoef = output.Fail;
    
    output.Data = struct();
    output.Data.Dec = sum(indices);
    output.Data.DecPrc = sum(indices) / params.Obs;
    output.Data.FrmAvgDec = returns_december_avg;
    output.Data.FrmAvgOth = returns_other_avg;
    output.Data.FrmSpr = spread;
    output.Data.H0Prc = percentile_h0;
    output.Data.H0Spr = spread_h0;

end

function res = test_digits_conformity(par,frm,ret,dec)

    a = par.A;

    ret = floor(abs(ret) .* (10 ^ dec));

    fst_ret = ret(ret >= 1);
    fst_n = numel(fst_ret);

    fst_the_dgts = (1:9).';
    fst_the_p = log10(1 + (1 ./ fst_the_dgts));

    fst_emp_dgts = (10 .^ ((floor(log10(fst_ret)) .* -1))) .* fst_ret;
    fst_emp_dgts = str2double(cellstr(num2str(fst_emp_dgts)));
    fst_emp_dgts = fst_emp_dgts - rem(fst_emp_dgts,1);
    fst_emp_p = histcounts(fst_emp_dgts,[fst_the_dgts; Inf]).' ./ fst_n;
    
    fst_chi2 = fst_n * sum(((fst_emp_p - fst_the_p) .^ 2) ./ fst_the_p);
    fst_pval = max([0 min([(1 - chi2cdf(fst_chi2,8)) 1])]);
    
    lst_ret = ret(ret >= 10);
    lst_n = numel(lst_ret);

    lst_the_dgts = (0:9).';
    lst_the_p = repmat(0.1,10,1);

    lst_emp_dgts = mod(lst_ret,10);
    lst_emp_p = histcounts(lst_emp_dgts,[lst_the_dgts; Inf]).' ./ lst_n;
    
    lst_chi2 = lst_n * sum(((lst_emp_p - lst_the_p) .^ 2) ./ lst_the_p);
    lst_pval = max([0 min([(1 - chi2cdf(lst_chi2,9)) 1])]);
    
    fail_fst = (fst_pval < a);
    fail_lst = (lst_pval < a);
    
    res = struct();
    res.Type = 'DC';
    res.Flags = 2;
    res.Par = par;
    res.Frm = frm;
    res.Fail = fail_fst | fail_lst;
    res.FailCoef = (fail_fst + fail_lst) / 2;

    res.Data = struct();
    res.Data.FstFail = fail_fst;
    res.Data.FstChi2 = fst_chi2;
    res.Data.FstEmpP = fst_emp_p;
    res.Data.FstTheP = fst_the_p;
    res.Data.FstPVal = fst_pval;
    res.Data.LstFail = fail_lst;
    res.Data.LstChi2 = lst_chi2;
    res.Data.LstEmpP = lst_emp_p;
    res.Data.LstTheP = lst_the_p;
    res.Data.LstPVal = lst_pval;

end

function output = test_discontinuity_at_zero(params,firm_name,returns)

    [bins,edges] = histcounts(returns,'BinWidth',0.005);
    edges_zero = find(edges == 0);

    x1 = bins(edges_zero - 2);
    p1 = x1 / params.Obs;

    x2 = bins(edges_zero - 1);
    p2 = x2 / params.Obs;

    x3 = bins(edges_zero);
    p3 = x3 / params.Obs;
    
    s = params.Obs * p2 * (1 - p2);
    
    if (s < 25)
        x2 = x2 - 0.5;
        p2 = x2 / params.Obs;
    end
    
    diff = x2 - mean([x1 x3]);
	diff_var = (params.Obs * p2 * (1 - p2)) + (0.25 * params.Obs * (p1 + p3) * (1 - p1 - p3)) + (params.Obs * p2 * (p1 + p3));

    z = diff / sqrt(diff_var);
    pval = 2 * normcdf(z);
    
    output = struct();

    output.Type = 'Discontinuity At Zero';
    output.Flags = 1;
    output.Par = params;
    output.Frm = firm_name;

    output.Fail = (diff < 0) & (pval < params.A);
    output.FailCoef = output.Fail;

    output.Data = struct();
    output.Data.Bins = bins;
    output.Data.Diff = diff;
    output.Data.DiffVar = diff_var;
    output.Data.Edges = edges;
    output.Data.IdxZ = edges_zero;
    output.Data.PVal = pval;
    output.Data.ZScore = z;

end

function [res,returns_fitted] = test_low_correlation(par,firm_name,returns,returns_h0,returns_other,style_factors,style_factors_comb,style_factors_names,r2_switch)

    a = par.A;
    n = par.Obs;
    t = par.Dates;

    sf_coms_k = size(style_factors_comb,2);
    sf_coms_n = size(style_factors_comb,1);

    o = ones(n,1);
    df = (n - 1) / (n - sf_coms_k - 1);
    returns_h0 = returns_h0(:,randperm(size(returns_h0,2),100));

    idx_ewi = mean(returns_other,2);
    idx_mdl = fitlm(idx_ewi,returns);
    idx_pval = idx_mdl.Coefficients{2,4};

    if (r2_switch)
        swi_fst = floor(0.1 * n);
        swi_lst = n - swi_fst;
        swi_pts = swi_fst:swi_lst;
        swi_pts_len = numel(swi_pts);

        swi_es2 = sum((returns - mean(returns)) .^ 2);
        swi_f = zeros(swi_pts_len,1);

        for i = 1:swi_pts_len
            swi_pt = swi_pts(i);

            swi_y1 = returns(1:swi_pt);
            swi_r1 = swi_y1 - mean(swi_y1);

            swi_y2 = returns((swi_pt + 1):n);
            swi_r2 = swi_y2 - mean(swi_y2);

            swi_us2 = sum(([swi_r1; swi_r2]) .^ 2);
            swi_f(i) = (swi_es2 - swi_us2) / (swi_us2 / (n - 2));
        end

        [~,swi_max_f] = max(swi_f);
        cp = swi_max_f + swi_fst;

        swi_iter = sf_coms_n * (sf_coms_n - 1);

        swi_cmbs = cell(swi_iter,2);
        swi_off = 1;

        for i = 1:sf_coms_n
            sf_cmb_i = style_factors_comb(i,:);

            for j = 1:sf_coms_n
                if (i == j)
                    continue;
                end

                sf_cmb_j = style_factors_comb(j,:);
                sf_curr = [style_factors{1:(cp-1),sf_cmb_i}; style_factors{cp:n,sf_cmb_j}];

                swi_cmbs(swi_off,:) = {sf_curr [style_factors_names(sf_cmb_i) {'>'} style_factors_names(sf_cmb_j)]};
                swi_off = swi_off + 1;
            end
        end

        swi_r2 = zeros(swi_iter,1);

        for i = 1:swi_iter
            [~,~,~,~,stat] = regress(returns,[o swi_cmbs{i,1}]);
            r2 = 1 - (1 - stat(1)) * df;

            swi_r2(i) = r2;
        end
        
        [r2_max,r2_max_idx] = max(swi_r2);
        max_nam = swi_cmbs{r2_max_idx,2};
        max_mod = fitlm(swi_cmbs{r2_max_idx,1},returns);
        
        ret_h0_max_r2s = -Inf(100,1);

        for i = 1:100
            ret_h0_curr = returns_h0(:,i);
            
            for j = 1:swi_iter
                [~,~,~,~,stat] = regress(ret_h0_curr,[o swi_cmbs{j,1}]);
                r2 = 1 - (1 - stat(1)) * df;

                if (r2 > ret_h0_max_r2s(i))
                    ret_h0_max_r2s(i) = r2;
                end
            end
        end

        ret_h0_prc = prctile(ret_h0_max_r2s,((1 - a) * 100));
        
        cp =  t(cp);
    else
        max_regs = cell(sf_coms_n,2);

        for i = 1:sf_coms_n
            sf_cmb = style_factors_comb(i,:);

            [~,~,~,~,stat] = regress(returns,[o style_factors(:,sf_cmb)]);
            r2 = 1 - (1 - stat(1)) * df;

            max_regs(i,:) = {r2 sf_cmb};
        end

        [r2_max,r2_max_idx] = max([max_regs{:,1}]);
        max_comb = max_regs{r2_max_idx,2};
        max_nam = style_factors_names(max_comb);
        max_mod = fitlm(style_factors(:,max_comb),returns);
        
        ret_h0_max_r2s = -Inf(100,1);

        for i = 1:100
            ret_h0_curr = returns_h0(:,i);
            
            for j = 1:sf_coms_n
                sf_cmb = style_factors_comb(j,:);

                [~,~,~,~,stat] = regress(ret_h0_curr,[o style_factors(:,sf_cmb)]);
                r2 = 1 - (1 - stat(1)) * df;

                if (r2 > ret_h0_max_r2s(i))
                    ret_h0_max_r2s(i) = r2;
                end
            end
        end

        ret_h0_prc = prctile(ret_h0_max_r2s,((1 - a) * 100));

        cp = [];
    end

    fail_idx = (idx_pval >= a);
    fail_max = (r2_max < ret_h0_prc);
    
    res = struct();
    res.Type = 'LC';
    res.Flags = 2;
    res.Par = par;
    res.Frm = firm_name;
    res.Fail = fail_idx | fail_max;
    res.FailCoef = (fail_idx + fail_max) / 2;

    res.Data = struct();
    res.Data.IdxFail = fail_idx;
    res.Data.IdxEWI = idx_ewi;
    res.Data.IdxMod = idx_mdl;
    res.Data.MaxFail = fail_max;
    res.Data.MaxComb = max_nam;
    res.Data.MaxCP = cp;
    res.Data.MaxH0 = ret_h0_max_r2s;
    res.Data.MaxMod = max_mod;
    res.Data.MaxTest = ret_h0_prc;
    res.Data.MaxVal = r2_max;
    
    returns_fitted = max_mod.Fitted;

end

function output = test_serial_correlation(params,firm_name,returns,returns_fitted)

    y = returns(2:end);
    x1 = returns(1:end-1);
    x2 = (1 - (returns_fitted(1:end-1) > mean(returns_fitted))) .* x1;
    
    mdl_conditional = fitlm([x1 x2],y);
    b_conditional = mdl_conditional.Coefficients{3,1};
    pval_conditional = mdl_conditional.Coefficients{3,4};

    mdl_unconditional = fitlm(x1,y);
    b_unconditional = mdl_unconditional.Coefficients{2,1};
    pval_unconditional = mdl_unconditional.Coefficients{2,4};

    failure_conditional = (b_conditional > 0) & (pval_conditional < params.A);
    failure_unconditional = (b_unconditional > 0) & (pval_unconditional < params.A);
    
    output = struct();

    output.Type = 'Serial Correlation';
    output.Flags = 2;
    output.Par = params;
    output.Frm = firm_name;

    output.Fail = failure_conditional | failure_unconditional;
    output.FailCoef = (failure_conditional + failure_unconditional) / 2;
    
    output.Data = struct();
    output.Data.ConFail = failure_conditional;
    output.Data.ConMod = mdl_conditional;
    output.Data.UncFail = failure_unconditional;
    output.Data.UncMod = mdl_unconditional;

end