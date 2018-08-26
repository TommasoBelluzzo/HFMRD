% [INPUT]
% file = A string representing the full path to the Excel spreadsheet containing the dataset.
%
% [OUTPUT]
% data = A structure containing the parsed dataset.

function data = parse_dataset(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('file',@(x)validateattributes(x,{'char','string'},{'scalartext','nonempty'}));
    end
    
    p.parse(varargin{:});
    res = p.Results;

    data = parse_dataset_internal(res.file);

end

function data = parse_dataset_internal(file)

    if (exist(file,'file') == 0)
        error('The dataset file does not exist.');
    end

    [file_stat,file_shts,file_fmt] = xlsfinfo(file);

    if (isempty(file_stat) || (ispc() && ~strcmp(file_fmt,'xlOpenXMLWorkbook')))
        error('The dataset file is not a valid Excel spreadsheet.');
    end

    shts_len = length(file_shts);
    
    if (shts_len < 2)
        error('The dataset does not contain all the required sheets.');
    end
    
    if (shts_len > 3)
        error('The dataset contains unnecessary sheets.');
    end
    
    file_shts = strtrim(file_shts);
    
    if (~isequal(file_shts(1:2),{'Returns' 'Style Factors'}))
        error('The dataset contains invalid (wrong name) or misplaced (wrong order) sheets.');
    end
    
    rets = parse_table(file,1,'Returns',true);

    if (width(rets) < 6)
        error('The ''Returns'' table must contain at least the following series: observations dates, benchmark returns, risk-free rates and the returns of three firms to analyze.');
    end

    if (~strcmp(rets.Properties.VariableNames(2),'BM'))
        error('The second column of the ''Returns'' table must be called ''BM'' and must contain the benchmark returns.');
    end

    if (~strcmp(rets.Properties.VariableNames(3),'RF'))
        error('The third column of the ''Returns'' table must be called ''RF'' and must contain the risk-free rates.');
    end
    
    n = height(rets);

    if (n < 120)
        error('The dataset must contain at least 120 observations in order to run consistent calculations.');
    end
    
    if (any(ismissing(rets)))
        error('The ''Returns'' table contains invalid or missing values.');
    end
    
    dates_str = cellstr(datestr(rets{:,1},'mm/yyyy'));
    dates_num = datenum(rets{:,1});
    dates_beg = dates_num(1);
    dates_end = dates_num(end);
    rets.Date = [];
    
    bm = rets{:,1};
    rf = rets{:,2};

    frms = numel(rets.Properties.VariableNames) - 2;
    frms_nam = rets.Properties.VariableNames(3:end);
    frms_ret = rets{:,3:end};
    
    sf = parse_table(file,2,'Style Factors',true);
    
    if (width(sf) < 3)
        error('The ''Style Factors'' table must contain at least 3 style factors.');
    end
    
    if (height(sf) ~= n)
        error('The number of observations in the ''Returns'' table and in the ''Style Factors'' table must be identical.');
    end
    
    if (any(strcmp(rets.Properties.VariableNames,'MRKEXC')))
        error('The ''Style Factors'' table contains a column called ''MRKEXC'', which is a reserved name for the market excess returns.');
    end

    if (any(ismissing(sf)))
        error('The ''Style Factors'' table contains invalid or missing values.');
    end
    
    if ((datenum(sf.Date(1)) ~=  dates_beg) || (datenum(sf.Date(end)) ~=  dates_end) || (size(sf,1) ~= n))
        error('The ''Returns'' table and the ''Style Factors'' table observation dates are mismatching.');
    end

    mrk_exc = bm - rf;
    
    sf.Date = [];
    sf = [table(mrk_exc,'VariableNames',{'MRKEXC'}) sf];
    
    if (shts_len == 3)
        if (~strcmp(file_shts(3),'Groups'))
            error('The dataset contains invalid (wrong name) or misplaced (wrong order) sheets.');
        end
        
        grps = parse_table(file,3,'Groups',false);

        if (~isequal(grps.Properties.VariableNames,rets.Properties.VariableNames(3:end)))
            error('The ''Returns'' table and the ''Groups'' table fniirms are mismatching.');
        end
        
        if (height(grps) ~= 1)
            error('The ''Groups'' table must contain one row of values.');
        end
        
        grps_vals = grps{:,:};
        
        if (any(ismissing(grps_vals)) || any(grps_vals < 1))
            error('The ''Groups'' table contains invalid or missing values.');
        end
        
        grps_max = max(grps_vals);
        grps_seq = 1:grps_max;
        grps_cnt = zeros(grps_max,1);
        
        for i = grps_seq
            grps_cnt(i) = numel(grps_vals(grps_vals == i));
        end
        
        grps_mis = (grps_cnt == 0);

        if (any(grps_mis))
            grps_mis = sprintfc(' %d',grps_seq(grps_mis));
            error('The following groups are not defined in the ''Groups'' table:%s.',strcat(grps_mis{:}));
        end
        
        if (any(grps_cnt < 3))
            error('Each group defined in the ''Groups'' table must contain at least 3 firms.');
        end
    else
        grps_vals = ones(1,frms);
    end
    
    data = struct();
    data.BM = bm;
    data.BMAnn = (prod(1 + bm) ^ (12 / n)) - 1;
    data.BMAvg = mean(bm);
    data.DatesNum = dates_num;
    data.DatesStr = dates_str;
    data.Frms = frms;
    data.FrmsNam = frms_nam;
    data.FrmsRet = frms_ret;
    data.Grps = grps_vals;
    data.MAR = mean(rf) + (2 * std(rf));
    data.ME = mrk_exc;
    data.MEAnn = (prod(1 + mrk_exc) ^ (12 / n)) - 1;
    data.MEAvg = mean(mrk_exc);
    data.RF = rf;
    data.RFAnn =  (prod(1 + rf) ^ (12 / n)) - 1;
    data.RFAvg =  mean(rf);
    data.Obs = n;
    data.PerAdj = sqrt(12);
    data.RiskAve = (log(mean(1 + bm)) - log(mean(1 + rf))) / var(1 + bm);
    data.SF = sf;

end

function res = parse_table(file,sht,name,ts)

    if (verLessThan('Matlab','9.1'))
        res = readtable(file,'Sheet',sht);
        
        if (~all(cellfun(@isempty,regexp(res.Properties.VariableNames,'^Var\d+$','once'))))
            error(['The ''' name ''' table contains unnamed columns.']);
        end
        
        if (ts)
            if (~strcmp(res.Properties.VariableNames(1),'Date'))
                error(['The first column of the ''' name ''' table must be called ''Date'' and must contain the time series dates.']);
            end

            res.Date = datetime(res.Date,'InputFormat','MM/yyyy');
            res_vars = varfun(@class,res,'OutputFormat','cell');

            if (~all(strcmp(res_vars(2:end),'double')))
                error(['The ''' name ''' table contains invalid or missing values.']);
            end
        else
            res_vars = varfun(@class,res,'OutputFormat','cell');
            
            if (~all(strcmp(res_vars,'double')))
                error(['The ''' name ''' table contains invalid or missing values.']);
            end
        end
    else
        opts = detectImportOptions(file,'Sheet',sht);
        
        if (~all(cellfun(@isempty,regexp(opts.VariableNames,'^Var\d+$','once'))))
            error(['The ''' name ''' table contains unnamed columns.']);
        end
        
        if (ts)
            if (~strcmp(opts.VariableNames(1),'Date'))
                error(['The first column of the ''' name ''' table must be called ''Date'' and must contain the time series dates.']);
            end

            opts = setvartype(opts,[{'datetime'} repmat({'double'},1,numel(opts.VariableNames)-1)]);
            opts = setvaropts(opts,'Date','InputFormat','MM/yyyy');
        else
            opts = setvartype(opts,repmat({'double'},1,numel(opts.VariableNames)));
        end
        
        res = readtable(file,opts);
    end

end
