% [INPUT]
% data = A structure representing the dataset.
% tab  = A boolean indicating whether the tests result is returned as a summary table or a cell array of structures with detailed information (optional, default=true).
% a    = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.10).
% sims = An integer representing the number of Monte Carlo simulations to perform (optional, default=10000).
% nsf  = An integer [1,3] representing the maximum number of style factors that can be simultaneously used in regression models (optional, default=3).
% swi  = A boolean indicating whether to use a switch of style factors for the calculation of the maximum R2 in the low correlation test (optional, default=false).
% dec  = An integer [2,6] representing the number of decimal places to consider when performing digit tests (optional, default=4).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
%
% [OUTPUT]
% td   = An instance containing the test results, whose type depends on the "tab" input parameter:
%         - true: a summary table that displays the outcome and the anomaly score of each test;
%         - false: a cell array of structures with the quantitative details of each test.

function td = execute_tests(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'struct'},{'nonempty'}));
        p.addOptional('tab',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
        p.addOptional('a',0.10,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('sims',10000,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1000}));
        p.addOptional('nsf',3,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1,'<=',6}));
        p.addOptional('swi',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
        p.addOptional('dec',4,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',2,'<=',6}));
    end
    
    p.parse(varargin{:});
    res = p.Results;
    data = res.data;
    tab = res.tab;
    a = res.a;
    sims = res.sims;
    nsf = res.nsf;
    swi = res.swi;
    dec = res.dec;

    if (nsf > width(data.SF))
        error('The maximum number of style factors must be less than or equal to the number of style factors plus one.');
    end

    td = execute_tests_internal(data,tab,a,sims,nsf,swi,dec);

end

function td = execute_tests_internal(data,tab,a,sims,nsf,swi,dec)

    par = struct();
    par.A = a;
    par.Dates = data.DatesNum;
    par.Obs = data.Obs;
    par.Sims = sims;
    
    sf = data.SF;
    sf_coms = nchoosek(1:width(sf),nsf);

    test_int = cell(7,data.Frms);
    
    for i = 1:data.Frms
        frm = data.FrmsNam{i};
        
        ret = data.FrmsRet(:,i);
        ret_h0 = round(normrnd(mean(ret),std(ret),[data.Obs sims]),4);
        ret_oth = data.FrmsRet(:,(~strcmp(data.FrmsNam,data.FrmsNam(i)) & (data.Grps == data.Grps(i))));

        [res_lco,ret_fit] = test_low_correlation(par,frm,ret,ret_h0,ret_oth,sf,sf_coms,swi);

        test_int{1,i} = res_lco;
        test_int{2,i} = test_serial_correlation(par,frm,ret,ret_fit);
        test_int{3,i} = test_bias_ratio(par,frm,ret,ret_h0);
        test_int{4,i} = test_december_spike(par,frm,ret,ret_h0);
        test_int{5,i} = test_discontinuity_at_zero(par,frm,ret);
        test_int{6,i} = test_digits_conformity(par,frm,ret,dec);
        test_int{7,i} = test_data_quality(par,frm,ret,ret_h0);
    end

    if (tab)
        td = cell(data.Frms,1);

        for j = 1:data.Frms
            fail = false(8,1);
            coef = zeros(8,1);

            for i = 1:7
                test_curr = test_int{i,j};
                fail(i) = test_curr.Fail;
                coef(i) = test_curr.FailCoef;
            end

            fail(8) = any(fail(1:7));
            coef(8) = sum(coef(1:7));

            td{j} = table(fail,coef,'VariableNames',{'Failure' 'Coefficient'});
        end

        td = table(td{:},'VariableNames',data.FrmsNam);
        td.Properties.RowNames = {'Low Correlation' 'Serial Correlation' 'Bias Ratio' 'December Spike' 'Discontinuity At Zero' 'Digits Conformity' 'Data Quality' 'Total'};
    else
        td = test_int;
    end

end

function res = test_bias_ratio(par,frm,ret,ret_h0)

    a = par.A;

    ub = std(ret);
    lb = -ub;
    br = sum((ret >= 0) & (ret <= ub)) / sum((ret >= lb) & (ret < 0));

    ret_h0_len = size(ret_h0,2);
    ret_h0_brs = zeros(ret_h0_len,1);
    
    for i = 1:ret_h0_len
        ret_h0_curr = ret_h0(:,i);
        
        ret_h0_ub = std(ret_h0_curr);
        ret_h0_lb = - ret_h0_ub;
        ret_h0_brs(i) = sum((ret_h0_curr >= 0) & (ret_h0_curr <= ret_h0_ub)) / sum((ret_h0_curr >= ret_h0_lb) & (ret_h0_curr < 0));
    end
    
    pval = sum(ret_h0_brs >= (br - 1e-8)) / ret_h0_len;
    pval = max([0 min([pval 1])]);
    
    res = struct();
    res.Type = 'BR';
    res.Flags = 1;
    res.Par = par;
    res.Frm = frm;
    res.Fail = (br >= 2.5) & (pval < a);
    res.FailCoef = res.Fail;
    
    res.Data = struct();
    res.Data.BR = br;
    res.Data.H0 = ret_h0_brs;
    res.Data.PVal = pval;

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

function res = test_december_spike(par,frm,ret,ret_h0)

    a = par.A;
    n = par.Obs;
    t = par.Dates;

    is_dec = (month(t) == 12);
    
    ret_dec = ret(is_dec,:);
    ret_dec_avg = mean(ret_dec);
    ret_oth = ret(~is_dec,:);
    ret_oth_avg = mean(ret_oth);
    ret_spr = ret_dec_avg - ret_oth_avg;
    
    ret_h0_dec = ret_h0(is_dec,:);
    ret_h0_dec_avg = mean(ret_h0_dec);
    ret_h0_oth = ret_h0(~is_dec,:);
    ret_h0_oth_avg = mean(ret_h0_oth);
    ret_h0_spr = ret_h0_dec_avg - ret_h0_oth_avg;
    ret_h0_prc = prctile(ret_h0_spr,((1 - a) * 100));
    
    res = struct();
    res.Type = 'DS';
    res.Flags = 1;
    res.Par = par;
    res.Frm = frm;
    res.Fail = (ret_spr > ret_h0_prc);
    res.FailCoef = res.Fail;
    
    res.Data = struct();
    res.Data.Dec = sum(is_dec);
    res.Data.DecPrc = sum(is_dec) / n;
    res.Data.FrmAvgDec = ret_dec_avg;
    res.Data.FrmAvgOth = ret_oth_avg;
    res.Data.FrmSpr = ret_spr;
    res.Data.H0Prc = ret_h0_prc;
    res.Data.H0Spr = ret_h0_spr;

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

function res = test_discontinuity_at_zero(par,frm,ret)

    a = par.A;
    n = par.Obs;

    [bins,edgs] = histcounts(ret,'BinWidth',0.005);
    edg_zero = find(edgs == 0);

    x1 = bins(edg_zero - 2);
    p1 = x1 / n;

    x2 = bins(edg_zero - 1);
    p2 = x2 / n;

    x3 = bins(edg_zero);
    p3 = x3 / n;
    
    s = n * p2 * (1 - p2);
    
    if (s < 25)
        x2 = x2 - 0.5;
        p2 = x2 / n;
    end
    
    d = x2 - mean([x1 x3]);
	d_var = (n * p2 * (1 - p2)) + (0.25 * n * (p1 + p3) * (1 - p1 - p3)) + (n * p2 * (p1 + p3));

    z = d / sqrt(d_var);
    pval = 2 * normcdf(z);
    
    res = struct();
    res.Type = 'DZ';
    res.Flags = 1;
    res.Par = par;
    res.Frm = frm;
    res.Fail = (d < 0) & (pval < a);
    res.FailCoef = res.Fail;

    res.Data = struct();
    res.Data.Bins = bins;
    res.Data.Diff = d;
    res.Data.DiffVar = d_var;
    res.Data.Edges = edgs;
    res.Data.IdxZ = edg_zero;
    res.Data.PVal = pval;
    res.Data.ZScore = z;

end

function [res,ret_fit] = test_low_correlation(par,frm,ret,ret_h0,ret_oth,sf,sf_cmbs,swi)

    a = par.A;
    n = par.Obs;
    t = par.Dates;

    sf_coms_k = size(sf_cmbs,2);
    sf_coms_n = size(sf_cmbs,1);

    o = ones(n,1);
    df = (n - 1) / (n - sf_coms_k - 1);
    ret_h0 = ret_h0(:,randperm(size(ret_h0,2),100));

    idx_ewi = mean(ret_oth,2);
    idx_mdl = fitlm(idx_ewi,ret);
    idx_pval = idx_mdl.Coefficients{2,4};

    if (swi)
        swi_fst = floor(0.1 * n);
        swi_lst = n - swi_fst;
        swi_pts = swi_fst:swi_lst;
        swi_pts_len = numel(swi_pts);

        swi_es2 = sum((ret - mean(ret)) .^ 2);
        swi_f = zeros(swi_pts_len,1);

        for i = 1:swi_pts_len
            swi_pt = swi_pts(i);

            swi_y1 = ret(1:swi_pt);
            swi_r1 = swi_y1 - mean(swi_y1);

            swi_y2 = ret((swi_pt + 1):n);
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
            sf_cmb_i = sf_cmbs(i,:);

            for j = 1:sf_coms_n
                if (i == j)
                    continue;
                end

                sf_cmb_j = sf_cmbs(j,:);
                sf_curr = [sf{1:(cp-1),sf_cmb_i}; sf{cp:n,sf_cmb_j}];

                swi_cmbs(swi_off,:) = {sf_curr [sf.Properties.VariableNames(sf_cmb_i) {'>'} sf.Properties.VariableNames(sf_cmb_j)]};
                swi_off = swi_off + 1;
            end
        end

        swi_r2 = zeros(swi_iter,1);

        for i = 1:swi_iter
            [~,~,~,~,stat] = regress(ret,[o swi_cmbs{i,1}]);
            r2 = 1 - (1 - stat(1)) * df;

            swi_r2(i) = r2;
        end
        
        [r2_max,r2_max_idx] = max(swi_r2);
        max_nam = swi_cmbs{r2_max_idx,2};
        max_mod = fitlm(swi_cmbs{r2_max_idx,1},ret);
        
        ret_h0_max_r2s = -Inf(100,1);

        for i = 1:100
            ret_h0_curr = ret_h0(:,i);
            
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
            sf_cmb = sf_cmbs(i,:);

            [~,~,~,~,stat] = regress(ret,[o sf{:,sf_cmb}]);
            r2 = 1 - (1 - stat(1)) * df;

            max_regs(i,:) = {r2 sf_cmb};
        end

        [r2_max,r2_max_idx] = max([max_regs{:,1}]);
        max_comb = max_regs{r2_max_idx,2};
        max_nam = sf.Properties.VariableNames(max_comb);
        max_mod = fitlm(sf{:,max_comb},ret);
        
        ret_h0_max_r2s = -Inf(100,1);

        for i = 1:100
            ret_h0_curr = ret_h0(:,i);
            
            for j = 1:sf_coms_n
                sf_cmb = sf_cmbs(j,:);

                [~,~,~,~,stat] = regress(ret_h0_curr,[o sf{:,sf_cmb}]);
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
    res.Frm = frm;
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
    
    ret_fit = max_mod.Fitted;

end

function res = test_serial_correlation(par,frm,ret,ret_fit)

    a = par.A;

    i = ret_fit(1:end-1) > mean(ret_fit);

    y = ret(2:end);
    x1 = ret(1:end-1);
    x2 = (1 - i) .* x1;
    
    con_mdl = fitlm([x1 x2],y);
    con_b = con_mdl.Coefficients{3,1};
    con_pval = con_mdl.Coefficients{3,4};

    unc_mdl = fitlm(x1,y);
    unc_b = unc_mdl.Coefficients{2,1};
    unc_pval = unc_mdl.Coefficients{2,4};

    fail_con = (con_b > 0) & (con_pval < a);
    fail_unc = (unc_b > 0) & (unc_pval < a);
    
    res = struct();
    res.Type = 'SC';
    res.Flags = 2;
    res.Par = par;
    res.Frm = frm;
    res.Fail = fail_con | fail_unc;
    res.FailCoef = (fail_con + fail_unc) / 2;
    
    res.Data = struct();
    res.Data.ConFail = fail_con;
    res.Data.ConMod = con_mdl;
    res.Data.UncFail = fail_unc;
    res.Data.UncMod = unc_mdl;

end