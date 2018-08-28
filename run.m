warning('off','all');

close('all');
clearvars();
clc();

[path,~,~] = fileparts(mfilename('fullpath'));

if (~endsWith(path,filesep()))
    path = [path filesep()];
end

paths_base = genpath(path);
addpath(paths_base);

data = parse_dataset(fullfile(path,'\Datasets\Example.xlsx'));
td = execute_tests(data,false);

plot_data(data,true);
plot_results(td);

rmpath(paths_base);