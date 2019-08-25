% [INPUT]
% file        = A string representing the full path to the Excel spreadsheet containing the dataset.
% date_format = A string representing the date format used in the Excel spreadsheet (optional, default=MM/yyyy).
%
% [OUTPUT]
% data        = A structure containing the parsed dataset.

function data = parse_dataset(varargin)

    persistent ip;

    if (isempty(ip))
        ip = inputParser();
        ip.addRequired('file',@(x)validateattributes(x,{'char'},{'nonempty','size',[1,NaN]}));
        ip.addOptional('date_format','MM/yyyy',@(x)validateattributes(x,{'char'},{'nonempty','size',[1,NaN]}));
    end

    ip.parse(varargin{:});
    ipr = ip.Results;

    data = parse_dataset_internal(ipr.file,ipr.date_format);

end

function data = parse_dataset_internal(file,date_format)

    try
        datetime('now','InputFormat',date_format);
    catch
        error('The specified date format is invalid.');
    end

    if (exist(file,'file') == 0)
        error('The dataset file does not exist.');
    end

    if (ispc())
        [file_status,file_sheets,file_format] = xlsfinfo(file);
        
        if (isempty(file_status) || ~strcmp(file_format,'xlOpenXMLWorkbook'))
            error('The dataset file is not a valid Excel spreadsheet.');
        end
    else
        [file_status,file_sheets] = xlsfinfo(file);
        
        if (isempty(file_status))
            error('The dataset file is not a valid Excel spreadsheet.');
        end
    end
    
    if (~strcmp(file_sheets{1},'Returns'))
        error('The first sheet of the dataset file must be the ''Returns'' one.');
    end
    
    if (~ismember('Style Factors',file_sheets))
        error('The dataset does not contain the required ''Style Factors'' sheet.');
    end
    
    tab_returns = parse_table(file,1,'Returns',date_format);

    if (width(tab_returns) < 6)
        error('The ''Returns'' table must contain at least the following series: observations dates, benchmark returns, risk-free rates and the returns of 3 firms to analyze.');
    end

    if (~strcmp(tab_returns.Properties.VariableNames(2),'BM'))
        error('The second column of the ''Returns'' table must be called ''BM'' and contain the benchmark returns.');
    end

    if (~strcmp(tab_returns.Properties.VariableNames(3),'RF'))
        error('The third column of the ''Returns'' table must be called ''RF'' and contain the risk-free rates.');
    end
    
    t = height(tab_returns);

    if (t < 126)
        error('The dataset must contain at least 126 observations (half of a business year) in order to run consistent calculations.');
    end
    
    if (any(any(ismissing(tab_returns))))
        error('The ''Returns'' table contains invalid or missing values.');
    end
    
    dates_str = cellstr(datestr(tab_returns{:,1},'mm/yyyy'));
    dates_num = datenum(tab_returns{:,1});
    tab_returns.Date = [];
    
    benchmark = tab_returns{:,1};
    risk_free = tab_returns{:,2};
    market_excess = benchmark - risk_free;

    firms = numel(tab_returns.Properties.VariableNames) - 2;
    firm_names = tab_returns.Properties.VariableNames(3:end);
    firm_returns = tab_returns{:,3:end};
    
    groups_indices = ones(1,firms);
    style_factors = table([]);
    style_factors_names = [];
    
    for tab = {'Style Factors' 'Groups'}
        
        tab_index = find(strcmp(file_sheets,tab),1);

        switch (char(tab))

            case 'Style Factors'

                tab_style_factors = parse_table(file,tab_index,'Style Factors',date_format);

                if (width(tab_style_factors) < 3)
                    error('The ''Style Factors'' table must contain at least 3 time series.');
                end

                if (any(any(ismissing(tab_style_factors))))
                    error('The ''Style Factors'' sheet contains invalid or missing values.');
                end

                if ((size(tab_style_factors,1) ~= t) || any(datenum(tab_style_factors.Date) ~= dates_num))
                    error('The observation dates in ''Returns'' and ''Style Factors'' sheets are mismatching.');
                end

                if (any(strcmp(tab_style_factors.Properties.VariableNames,'MRKEXC')))
                    error('The ''Style Factors'' table contains a column called ''MRKEXC'', which is a reserved name for the market excess returns.');
                end

                tab_style_factors.Date = [];
                tab_style_factors = [table(market_excess,'VariableNames',{'MRKEXC'}) tab_style_factors]; %#ok<AGROW>

                style_factors = tab_style_factors{:,:};
                style_factors_names = tab_style_factors.Properties.VariableNames;

            case 'Groups'

                if (~isempty(tab_index))
                    groups = parse_table(file,tab_index,'Groups',date_format);

                    if (~isequal(groups.Properties.VariableNames,tab_returns.Properties.VariableNames(3:end)))
                        error('The firms defined in ''Returns'' table and ''Groups'' table are mismatching.');
                    end

                    if (height(groups) ~= 1)
                        error('The ''Groups'' table must contain only one row of values.');
                    end

                    groups_indices = groups{:,:};

                    if (any(ismissing(groups_indices)) || any(groups_indices < 1) || any(round(groups_indices) ~= groups_indices))
                        error('The ''Groups'' table contains invalid or missing values.');
                    end

                    groups_max = max(groups_indices);
                    groups_sequence = 1:groups_max;
                    groups_count = zeros(groups_max,1);

                    for i = groups_sequence
                        groups_count(i) = numel(groups_indices(groups_indices == i));
                    end

                    groups_undefined = (groups_count == 0);

                    if (any(groups_undefined))
                        groups_undefined = sprintfc(' %d',groups_sequence(groups_undefined));
                        error('The following groups are not defined in the ''Groups'' table:%s.',strcat(groups_undefined{:}));
                    end

                    if (numel(unique(groups_indices)) > 10)
                        error('A maximum of 10 groups can be defined the ''Groups'' table.');
                    end

                    if (any(groups_count < 3))
                        error('Each group defined in the ''Groups'' table must contain at least 3 firms.');
                    end
                end

        end
    end
    
    data = struct();
    
    data.Obs = t;
    data.Frms = firms;
    
    data.DatesNum = dates_num;
    data.DatesStr = dates_str;

    data.FirmNames = firm_names;
    data.FrmsRet = firm_returns;
    data.StyleFactors = style_factors;
    data.StyleFactorsNames = style_factors_names;
    data.Grps = groups_indices;

    data.Benchmark = benchmark;
    data.BenchmarkAnnualized = (prod(1 + benchmark) ^ (12 / t)) - 1;
    data.BenchmarkAverage = mean(benchmark);
    data.MarketExcess = market_excess;
    data.MarketExcessAnnualized = (prod(1 + market_excess) ^ (12 / t)) - 1;
    data.MarketExcessAverage = mean(market_excess);
    data.RiskFree = risk_free;
    data.RiskFreeAnnualized = (prod(1 + risk_free) ^ (12 / t)) - 1;
    data.RiskFreeAverage = mean(risk_free);

    data.MAR = mean(risk_free) + (2 * std(risk_free));
    data.PerAdj = sqrt(12);
    data.RiskAversion = (log(mean(1 + benchmark)) - log(mean(1 + risk_free))) / var(1 + benchmark);

end

function output = parse_table(file,sheet,name,date_format)

    if (verLessThan('Matlab','9.1'))
        output = readtable(file,'Sheet',sheet);
        
        if (~all(cellfun(@isempty,regexp(output.Properties.VariableNames,'^Var\d+$','once'))))
            error(['The ''' name ''' table contains unnamed columns.']);
        end

        if (strcmp(name,'Groups'))
            output_vars = varfun(@class,output,'OutputFormat','cell');
            
            if (~all(strcmp(output_vars,'double')))
                error(['The ''' name ''' table contains invalid or missing values.']);
            end
        else
            if (~strcmp(output.Properties.VariableNames(1),'Date'))
                error(['The first column of the ''' name ''' table must be called ''Date'' and must contain the observation dates.']);
            end
            
            output.Date = datetime(output.Date,'InputFormat',date_format);
            output_vars = varfun(@class,output,'OutputFormat','cell');
            
            if (~all(strcmp(output_vars(2:end),'double')))
                error(['The ''' name ''' table contains invalid or missing values.']);
            end
        end
    else
        options = detectImportOptions(file,'Sheet',sheet);
        
        if (~all(cellfun(@isempty,regexp(options.VariableNames,'^Var\d+$','once'))))
            error(['The ''' name ''' table contains unnamed columns.']);
        end

        if (strcmp(name,'Groups'))
            options = setvartype(options,repmat({'double'},1,numel(options.VariableNames)));
        else
            if (~strcmp(options.VariableNames(1),'Date'))
                error(['The first column of the ''' name ''' table must be called ''Date'' and must contain the observation dates.']);
            end

            options = setvartype(options,[{'datetime'} repmat({'double'},1,numel(options.VariableNames)-1)]);
            options = setvaropts(options,'Date','InputFormat',date_format);
        end

        output = readtable(file,options);
    end

end
