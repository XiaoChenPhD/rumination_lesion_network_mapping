% To extract the regional metrics from the rois 
% compare them in the Suzhou rumination dataset.
% 
% Xiao Chen 260508
%
% stable version
% Xiao Chen 260630
% chenxiaophd@gmail.com

%% initialization

clear; clc;

work_dir = 'xxxxxxxxxx'; %modify here
if ~exist(work_dir, "dir")
    mkdir(work_dir)
end
dpabiSurfTemplateDir = '/dpabi_dir/DPABISurf/SurfTemplates';
roi_dir = '~/rumination_lesion_network_mapping/maps/thresholded';
TopOutputDir = [work_dir, '/roi_signal_extracted'];
if ~exist(TopOutputDir, 'dir'); mkdir(TopOutputDir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Beijing data
% behavior_data = readtable(['xxxxxxxx/demographic_info_Beijing_v1.xlsx'], ...
%                             'ReadVariableNames', true);
% behavior_data = behavior_data(behavior_data.Handedness_label ~= 2, :);%remove left handedness
% sub_list = behavior_data.Serial_number;
% scale_score = behavior_data;
% scale_name_set = scale_score.Properties.VariableNames(6:end);
% 
% ConditionSet = {'IPCAS_rum', 'IPCAS_dis', 'PKUGE_rum', 'PKUGE_dis', ...
%               'PKUSIEMENS_rum', 'PKUSIEMENS_dis'}; % Beijing
% 
% condition_label = {'IPCAS_rum', 'IPCAS_dis', 'PKUGE_rum', 'PKUGE_dis', ...
%                    'PKUSIEMENS_rum', 'PKUSIEMENS_dis', ...
%                    'IPCAS_delta', 'IPCAS_relative_delta', ...
%                    'PKUGE_delta', 'PKUGE_relative_delta', ...
%                    'PKUSIEMENS_delta', 'PKUGE_relative_delta'}; %i

% Suzhou data
behavior_data = readtable('xxxxxxxx/sample_cleaned_merged_v5.csv', ...
                            'ReadVariableNames', true);
sub_list = behavior_data.Serial_Number;
Dx = behavior_data.Dx_selected_data;
age = behavior_data.Age;
sex = behavior_data.Sex;
sex = sex -1;
edu = behavior_data.Edu;
HeadMotion = [behavior_data.head_motion_rum, behavior_data.head_motion_dis, ...
              behavior_data.head_motion_rest];
scale_data_full = readtable('xxxxxxxx/demographic_info_v10.xlsx');
scale_score = scale_data_full(ismember(scale_data_full.Serial_Number, sub_list), :);
scale_name_set = scale_score.Properties.VariableNames(11:end);
vars = varfun(@isnumeric, scale_score, 'OutputFormat', 'uniform');
scale_score{:, vars}(scale_score{:, vars} == 999) = NaN;
ConditionSet = {'MDD_rum', 'MDD_dis', 'MDD_rest', 'HC_rum', 'HC_dis', 'HC_rest'};% Suzhou
condition_label = {'MDD_Rum', 'MDD_Dis', 'MDD_rest', 'HC_Rum', 'HC_Dis', 'HC_rest', ...
                   'MDD_delta', 'HC_delta', ...
                   'MDD_relative_delta', 'HC_relative_delta'}; %i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HemisphereNameSet = {'lh', 'rh'};
metric_set = {'ALFF', 'ALFF', 'fALFF', 'fALFF', ...
              'ReHo', 'ReHo', 'DegreeCentrality', 'DegreeCentrality'};
affix_set = {'_FunSurfWC', '_FunSurfWglobalC', ...
                        '_FunSurfWC', '_FunSurfWglobalC', ...
                        '_FunSurfWCF', '_FunSurfWglobalCF', ...
                        '_FunSurfWCF', '_FunSurfWglobalCF'};
                    
roi_set_lh = {};
roi_set_rh = {'/positive_gt3p29_cluster1_n69_mask.func.gii', ...
              '/positive_gt3p29_cluster2_n50_mask.func.gii', ...
              '/positive_gt3p29_cluster3_n33_mask.func.gii'};
roi_label = {'sig_cluster1', 'sig_cluster2', 'sig_cluster3'};
MRI_label = roi_label; %j

%% extract regional metrics
IsMultipleLabel = 0;
GHeader = [];
CUTNUMBER = 1;

% left hemisphere
ROIDef = cell(length(roi_set_lh));
for iROI = 1:length(roi_set_lh)
    ROIDef{iROI} = [roi_dir, roi_set_lh{iROI}];
end
for iCondition = 1:length(ConditionSet)
    for i_metric = 1:length(metric_set)
        AllVolume = [work_dir, '/metrics_organized/', ConditionSet{iCondition}, '/', ...
            metric_set{i_metric}, affix_set{i_metric}, '_lh'];
        OutputName = [TopOutputDir, '/', ConditionSet{iCondition}, '_', ...
            metric_set{i_metric}, affix_set{i_metric},'_lh'];
        AMaskFilename = [dpabiSurfTemplateDir, '/fsaverage5_lh_cortex.label.gii'];
        [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, ...
                            AMaskFilename, IsMultipleLabel, GHeader, CUTNUMBER);
    end
end

% right hemisphere
ROIDef = cell(length(roi_set_rh));
for iROI = 1:length(roi_set_rh)
    ROIDef{iROI} = [roi_dir, roi_set_rh{iROI}];
end
for iCondition = 1:length(ConditionSet)
    for i_metric = 1:length(metric_set)
        AllVolume = [work_dir, '/metrics_organized/', ConditionSet{iCondition}, '/', ...
            metric_set{i_metric}, affix_set{i_metric}, '_rh'];
        OutputName = [TopOutputDir, '/', ConditionSet{iCondition}, '_', ...
            metric_set{i_metric}, affix_set{i_metric},'_rh'];
        AMaskFilename = [dpabiSurfTemplateDir, '/fsaverage5_rh_cortex.label.gii'];
        [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, ...
                            AMaskFilename, IsMultipleLabel, GHeader, CUTNUMBER);
    end
end

%% organize
% 3 rois: sig_cluster1, sig_cluster2, sig_cluster3
% Suzhou
% for each metrics, organize into a 1 x 6 cell:6 conditions
% 'MDD_rum', 'MDD_dis', 'MDD_rest', 'HC_rum', 'HC_dis', 'HC_rest'
% Beijing
% 'IPCAS_rum', 'IPCAS_dis', 'PKUGE_rum', 'PKUGE_dis', 'PKUSIEMENS_rum', 'PKUSIEMENS_dis'
% each cell has a n x 3 matrix: 3 rois
if ~exist([work_dir, '/roi_signal_organized/'], 'dir')
    mkdir([work_dir, '/roi_signal_organized/']);
end
for i_metric = 1:length(metric_set)
    roi_cells = cell(1, length(ConditionSet));
    for iCondition = 1:length(ConditionSet)
        ROISignals_rh = load([TopOutputDir, '/ROISignals_', ConditionSet{iCondition}, '_', ...
                    metric_set{i_metric}, affix_set{i_metric},'_rh.mat']);
        roi_matrix = ROISignals_rh.ROISignals;
        roi_cells{iCondition} = roi_matrix;
    end
    save([work_dir, '/roi_signal_organized/', ...
        metric_set{i_metric}, affix_set{i_metric}, '.mat'], 'roi_cells');
end

%% export .csv files for stats and plotting
% suzhou
for i_metric = 1:length(metric_set)
    load([work_dir, '/roi_signal_organized/', ...
        metric_set{i_metric}, affix_set{i_metric}, '.mat']);
    for i_cluster = 1:length(roi_set_rh)
        metric = [];
        condition = [];
        group = [];
        metric = [
            roi_cells{1}(:,i_cluster);
            roi_cells{2}(:,i_cluster);
            roi_cells{4}(:,i_cluster);
            roi_cells{5}(:,i_cluster)
        ];
        
        condition = [
            ones(size(roi_cells{1}, 1),1);
            2*ones(size(roi_cells{2}, 1),1);
            ones(size(roi_cells{4}, 1),1);
            2*ones(size(roi_cells{5}, 1),1)
        ];
        
        group = [
            ones(size(roi_cells{1}, 1) + size(roi_cells{2}, 1),1);
            2*ones(size(roi_cells{4}, 1) + size(roi_cells{5}, 1),1)
        ];
        
        T = table(metric, condition, group);
        writetable(T, [work_dir, '/roi_signal_organized/', ...
        metric_set{i_metric}, affix_set{i_metric}, '_cluster', num2str(i_cluster), '.csv']);
    end
end

%% beijing
for i_metric = 1:length(metric_set)
    load([work_dir, '/roi_signal_organized/', ...
        metric_set{i_metric}, affix_set{i_metric}, '.mat']);
    for i_cluster = 1:length(roi_set_rh)
        metric = [];
        condition = [];
        group = [];
        metric = [
            roi_cells{1}(:,i_cluster);
            roi_cells{2}(:,i_cluster);
            roi_cells{3}(:,i_cluster);
            roi_cells{4}(:,i_cluster);
            roi_cells{5}(:,i_cluster);
            roi_cells{6}(:,i_cluster)
        ];
        
        condition = [
            ones(size(roi_cells{1}, 1),1);
            2*ones(size(roi_cells{2}, 1),1);
            ones(size(roi_cells{3}, 1),1);
            2*ones(size(roi_cells{4}, 1),1);
            ones(size(roi_cells{5}, 1),1);
            2*ones(size(roi_cells{6}, 1),1);
        ];
        
        group = [
            ones(size(roi_cells{1}, 1) + size(roi_cells{2}, 1),1);
            2*ones(size(roi_cells{3}, 1) + size(roi_cells{4}, 1),1);
            3*ones(size(roi_cells{5}, 1) + size(roi_cells{6}, 1),1)
        ];
        
        T = table(metric, condition, group);
        writetable(T, [work_dir, '/roi_signal_organized/', ...
        metric_set{i_metric}, affix_set{i_metric}, '_cluster', num2str(i_cluster), '.csv']);
    end
end

%% construct MRI Phenotype
% for each regional metrics, construct a cell to store all the values
% Suzhou: 
% a 1 x 10 cell, each cell contains a n x 3 matrix:
% sig_cluster1, sig_cluster2, sig_cluster3
% cell 1-6: MDD_Rum, MDD_Dis, MDD_rest, HC_Rum, HC_Dis, HC_rest
% cell 5-6: MDD_Rum-MDD_Dis, HC_Rum-HC_Dis
% cell 7-8: (MDD_Rum-MDD_Dis)/MDD_Dis, (HC_Rum-HC_Dis)/HC_Dis
% Beijing
% a 1 x 12 cell
% cell 1-6: IPCAS_rum, IPCAS_dis, PKUGE_rum, PKUGE_dis, PKUSIEMENS_rum,
% PKUSIEMENS_dis
% cell 7-12: IPCAS_delta, IPCAS_relative_delta, PKUGE_delta,
% PKUGE_relative_delta, PKUSIEMENS_delta, PKUSIEMENS_relative_delta

if ~exist([work_dir, '/stats_rois/'], 'dir')
    mkdir([work_dir, '/stats_rois/']);
end


for i_metric = 1:length(metric_set)
    load([work_dir, '/roi_signal_organized/', ...
                metric_set{i_metric}, affix_set{i_metric}, '.mat']);
    MRI_Phenotype = roi_cells;

    % % Suzhou
    % MRI_Phenotype{7} = MRI_Phenotype{1} - MRI_Phenotype{2};
    % MRI_Phenotype{8} = MRI_Phenotype{4} - MRI_Phenotype{5};
    % MRI_Phenotype{9} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
    % MRI_Phenotype{10} = (MRI_Phenotype{4} - MRI_Phenotype{5})./MRI_Phenotype{5};
    % 
    % % do scale corr analyses, explorative
    % % Suzhou
    % % MDD patients
    % stats_corr = [];
    % temp_MDD = [];
    % for i = [1:3,7,9] %select MDD
    %     count = 0;
    %     for j = 1:size(MRI_Phenotype{1},2)
    %         for k = 11:size(scale_score,2) % from columne 11 to end
    %             y = MRI_Phenotype{i}(:,j);
    %             x = scale_score{Dx == 1, k};
    %             [r, p] = corr(x,y, 'Rows', 'complete');
    %             if p < 0.05
    %                 count = count + 1;
    %                 eval(['stats_corr.', condition_label{i}, ...
    %                     '{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k-10}]']);
    %                 eval(['stats_corr.', condition_label{i}, '{count, 2} = r']);
    % 
    %                 eval(['stats_corr.', condition_label{i}, '{count, 3} = p']);
    %                 %store the significant condition ids
    %                 temp_MDD(end+1,1) = i;
    %                 temp_MDD(end,2) = j;
    %                 temp_MDD(end,3) = k;
    %              end
    %         end
    %     end
    % end
    % 
    % % HC
    % temp_HC = [];
    % for i = [4:6,8,10] % select HC
    %     count = 0;
    %     for j = 1:size(MRI_Phenotype{1},2)
    %         for k = 11:size(scale_score,2) % from columne 11 to end
    % 
    %             y = MRI_Phenotype{i}(:,j);
    %             x = scale_score{Dx == 2, k};
    %             [r, p] = corr(x,y, 'Rows', 'complete');
    %             if p < 0.05
    %                 count = count + 1;
    %                 eval(['stats_corr.', condition_label{i}, ...
    %                     '{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k-10}]']);
    %                 eval(['stats_corr.', condition_label{i}, '{count, 2} = r']);
    %                 eval(['stats_corr.', condition_label{i}, '{count, 3} = p']);
    %                 %store the significant condition ids
    %                 temp_HC(end+1,1) = i;
    %                 temp_HC(end,2) = j;
    %                 temp_HC(end,3) = k;
    %              end
    %         end
    %     end
    % end
    
    % Beijing
    MRI_Phenotype{7} = MRI_Phenotype{1} - MRI_Phenotype{2};
    MRI_Phenotype{8} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
    MRI_Phenotype{9} = MRI_Phenotype{3} - MRI_Phenotype{4};
    MRI_Phenotype{10} = (MRI_Phenotype{3} - MRI_Phenotype{4})./MRI_Phenotype{4};
    MRI_Phenotype{11} = MRI_Phenotype{5} - MRI_Phenotype{6};
    MRI_Phenotype{12} = (MRI_Phenotype{5} - MRI_Phenotype{6})./MRI_Phenotype{6};


    % beijing
    stats_corr = [];
    temp_MDD = [];
    scatter_data = struct();
    plot_count = 0;
    for i = 1:size(MRI_Phenotype,2)
        count = 0;
        for j = 1:size(MRI_Phenotype{1},2)
            for k = 6:size(scale_score,2) % from columne 6 to end
                y = MRI_Phenotype{i}(:,j);
                x = scale_score{:, k};
                [r, p] = corr(x,y, 'Rows', 'complete');
                if p < 0.05
                    count = count + 1;
                    eval(['stats_corr.',condition_label{i},'{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k-5}]']);
                    eval(['stats_corr.',condition_label{i},'{count, 2} = r']);
                    eval(['stats_corr.',condition_label{i},'{count, 3} = p']);
                    %store the significant condition ids
                    temp_MDD(end+1,1) = i;
                    temp_MDD(end,2) = j;
                    temp_MDD(end,3) = k;
                end
            end
        end
    end

    save([work_dir, '/stats_rois/', ...
        metric_set{i_metric}, affix_set{i_metric}, '_corr.mat'], "stats_corr");
end

%% save out interested results, Suzhou
% organize
scale_interested = {'Rumination', 'Brooding', 'Reflection'};
metric_interested = {'sig_cluster1', 'sig_cluster2', 'sig_cluster3'};

i_metric = 6; % change here to get different metrics!!!!
% metric_set = {'ALFF', 'ALFF', 'fALFF', 'fALFF', ...
%               'ReHo', 'ReHo', 'DegreeCentrality', 'DegreeCentrality'};
load([work_dir, '/roi_signal_organized/', ...
            metric_set{i_metric}, affix_set{i_metric}, '.mat']);
MRI_Phenotype = roi_cells;
MRI_label = roi_label; %j

condition_label = {'MDD_Rum', 'MDD_Dis', 'MDD_rest', 'HC_Rum', 'HC_Dis', 'HC_rest', ...
                   'MDD_delta', 'HC_delta', ...
                   'MDD_relative_delta', 'HC_relative_delta'}; %i

MRI_Phenotype{7} = MRI_Phenotype{1} - MRI_Phenotype{2};
MRI_Phenotype{8} = MRI_Phenotype{4} - MRI_Phenotype{5};
MRI_Phenotype{9} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
MRI_Phenotype{10} = (MRI_Phenotype{4} - MRI_Phenotype{5})./MRI_Phenotype{5};

scatter_data_MDD = struct();
plot_count = 0;
for i = [1:3,7,9] %select MDD
    count = 0;
    for j = 1:size(MRI_Phenotype{1},2)
        for k = 11:size(scale_score,2) % from columne 11 to end
            y = MRI_Phenotype{i}(:,j);
            x = scale_score{Dx == 1, k};
            [r, p] = corr(x,y, 'Rows', 'complete');

            % save out the x and y for plots
            if ismember(MRI_label{j}, metric_interested) && ...
                    ismember(scale_name_set{k-10}, scale_interested)
                plot_count = plot_count + 1;
                scatter_data_MDD(plot_count).condition = condition_label{i};
                scatter_data_MDD(plot_count).MRI = MRI_label{j};
                scatter_data_MDD(plot_count).Scale = scale_name_set{k-10};
                scatter_data_MDD(plot_count).x = x;
                scatter_data_MDD(plot_count).y = y;
                scatter_data_MDD(plot_count).r = r;
                scatter_data_MDD(plot_count).p = p;
            end
        end
    end
end

% HC
scatter_data_HC = struct();
plot_count = 0;
for i = [4:6,8,10] % select HC
    count = 0;
    for j = 1:size(MRI_Phenotype{1},2)
        for k = 11:size(scale_score,2) % from columne 11 to end

            y = MRI_Phenotype{i}(:,j);
            x = scale_score{Dx == 2, k};
            [r, p] = corr(x,y, 'Rows', 'complete');
            % save out the x and y for plots
            if ismember(MRI_label{j}, metric_interested) && ...
                    ismember(scale_name_set{k-10}, scale_interested)
                plot_count = plot_count + 1;
                scatter_data_HC(plot_count).condition = condition_label{i};
                scatter_data_HC(plot_count).MRI = MRI_label{j};
                scatter_data_HC(plot_count).Scale = scale_name_set{k-10};
                scatter_data_HC(plot_count).x = x;
                scatter_data_HC(plot_count).y = y;
                scatter_data_HC(plot_count).r = r;
                scatter_data_HC(plot_count).p = p;
            end
        end
    end
end
save([work_dir, '/stats_rois/organized_', ...
    metric_set{i_metric}, affix_set{i_metric},'.mat'], ...
    "scatter_data_MDD", "scatter_data_HC");

%% save out file 

scatter_mat_file = [work_dir, '/stats_rois/organized_', ...
                    metric_set{i_metric}, affix_set{i_metric},'.mat'];
scatter_out_dir = [work_dir, '/stats_rois/txt_', ...
                      metric_set{i_metric}, affix_set{i_metric}];
if ~exist(scatter_out_dir, "dir")
    mkdir(scatter_out_dir);
end
load(scatter_mat_file, 'scatter_data_MDD', 'scatter_data_HC');

% MDD + HC together
n_plot = length(scatter_data_MDD);

for n = 1:9
    %------------------------
    % MDD
    %------------------------
    x1 = scatter_data_MDD(n).x;
    y1 = scatter_data_MDD(n).y;

    idx1 = ~isnan(x1) & ~isnan(y1);
    data1 = [x1(idx1), y1(idx1), ones(sum(idx1),1)];

    %------------------------
    % HC
    %------------------------
    x2 = scatter_data_HC(n).x;
    y2 = scatter_data_HC(n).y;

    idx2 = ~isnan(x2) & ~isnan(y2);
    data2 = [x2(idx2), y2(idx2), 2*ones(sum(idx2),1)];

    %------------------------
    % Merge
    %------------------------
    out_data = [data1; data2];

    % filename
    file_name = sprintf('%s_%s_%s.txt',...
        scatter_data_MDD(n).condition,...
        scatter_data_MDD(n).MRI,...
        scatter_data_MDD(n).Scale);

    writematrix(out_data, fullfile(scatter_out_dir,file_name),'Delimiter','tab');
end

%% save out interested results 
% Beijing
% organize
scale_interested = {'Rumination_total', 'Brooding', 'Reflection'};
metric_interested = {'sig_cluster1', 'sig_cluster2', 'sig_cluster3'};

i_metric = 6; % change here to get different metrics!!!!
% metric_set = {'ALFF', 'ALFF', 'fALFF', 'fALFF', ...
%               'ReHo', 'ReHo', 'DegreeCentrality', 'DegreeCentrality'};
load([work_dir, '/roi_signal_organized/', ...
            metric_set{i_metric}, affix_set{i_metric}, '.mat']);
MRI_Phenotype = roi_cells;

% Beijing
MRI_Phenotype{7} = MRI_Phenotype{1} - MRI_Phenotype{2};
MRI_Phenotype{8} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
MRI_Phenotype{9} = MRI_Phenotype{3} - MRI_Phenotype{4};
MRI_Phenotype{10} = (MRI_Phenotype{3} - MRI_Phenotype{4})./MRI_Phenotype{4};
MRI_Phenotype{11} = MRI_Phenotype{5} - MRI_Phenotype{6};
MRI_Phenotype{12} = (MRI_Phenotype{5} - MRI_Phenotype{6})./MRI_Phenotype{6};   
% beijing
stats_corr = [];
temp_MDD = [];
scatter_data = struct();
plot_count = 0;
for i = 1:size(MRI_Phenotype,2)
    count = 0;
    for j = 1:size(MRI_Phenotype{1},2)
        for k = 6:size(scale_score,2) % from columne 6 to end
            y = MRI_Phenotype{i}(:,j);
            x = scale_score{:, k};
            [r, p] = corr(x,y, 'Rows', 'complete');

            % save out the x and y for plots
            if ismember(MRI_label{j}, metric_interested) && ...
                    ismember(scale_name_set{k-5}, scale_interested)
                plot_count = plot_count + 1;
                scatter_data(plot_count).condition = condition_label{i};
                scatter_data(plot_count).MRI = MRI_label{j};
                scatter_data(plot_count).Scale = scale_name_set{k-5};
                scatter_data(plot_count).x = x;
                scatter_data(plot_count).y = y;
                scatter_data(plot_count).r = r;
                scatter_data(plot_count).p = p;
            end
        end
    end
end
save([work_dir, '/stats_rois/organized_', ...
    metric_set{i_metric}, affix_set{i_metric},'.mat'], "scatter_data");
%% write out the .txt files for plot

scatter_mat_file = [work_dir, '/stats_rois/organized_', ...
                    metric_set{i_metric}, affix_set{i_metric},'.mat'];
scatter_out_dir = [work_dir, '/stats_rois/txt_', ...
                      metric_set{i_metric}, affix_set{i_metric}];
if ~exist(scatter_out_dir, "dir")
    mkdir(scatter_out_dir);
end

loaded_scatter_data = load(scatter_mat_file, 'scatter_data');
scatter_data_to_write = loaded_scatter_data.scatter_data;

for i_row = 1:numel(scatter_data_to_write)
    condition = string(scatter_data_to_write(i_row).condition);

    if endsWith(condition, "_rum")
        x = scatter_data_to_write(i_row).x(:);
        y = scatter_data_to_write(i_row).y(:);
        out_data = [x, y];

        file_name = sprintf('%s_%s_%s.txt', ...
            char(condition), ...
            char(string(scatter_data_to_write(i_row).MRI)), ...
            char(string(scatter_data_to_write(i_row).Scale)));

        writematrix(out_data, fullfile(scatter_out_dir, file_name), ...
            'Delimiter', 'tab');
    end
end
