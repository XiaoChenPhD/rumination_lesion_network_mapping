% This script is to do the permutaion correction
% first get the real t maps, and calculate the z values
%
% Xiao Chen 251115
% Warpped up by Xiao Chen 261119
% chenxiaophd@gmail.com

%% initialization
clear; clc;

work_dir = 'xxxxxxxxxxxxx';
output_dir = [work_dir, '/z_results'];
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

group_set = {'MDD', 'HC'};
GSR_set = {'','_GSR'};

[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));
ProgramPath=fullfile(DPABIPath, 'DPARSF');
Cfg.MaskFileSurfLH = [DPABIPath, '/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii'];
Cfg.MaskFileSurfRH = [DPABIPath, '/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii'];

% tumor dataset
demog_info = readtable([work_dir, '/scale_data.xlsx'], 'sheet', 1);
subj_tumor_list = demog_info{:,1};
subj_tumor_list(1,:) =[]; %remove the example row
subj_tumor_list = cellfun(@lower, subj_tumor_list, 'UniformOutput', false); %lower_case

% questionnaires
questionnaire_table = demog_info{:, 90:95};
questionnaire_table(1,:) = [];

% age, sex, edu
other_covs = demog_info{:, 4:6};
other_covs(1,:) = [];

% add two control variables
questionnaire_table = [questionnaire_table, other_covs(:,1)]; % age
questionnaire_table = [questionnaire_table, other_covs(:,3)]; % edu
questionnaire_table = [questionnaire_table, randi([0,60],1,size(questionnaire_table,1))'];

questionnaire_set = {'reflection', 'brooding', 'rumination', 'BDI', ...
                                             'anxiety_s', 'anxiety_t', ...
                                             'age', 'edu', 'rand_vec'};
                                     
 %%
 TF_Flag = 'T';
 for i_questionnaire = 13
    SeedSeries = questionnaire_table(:, i_questionnaire);
    
    % set other covariates
    % tumor dataset
    if i_questionnaire < 4
        OtherCovariates = other_covs;
        %OtherCovariates = [other_covs, questionnaire_table(:, 4)];
    elseif i_questionnaire < 7
        OtherCovariates = other_covs;
    else
        OtherCovariates = []; % control variables
    end
    OtherCovariatesMatrix = OtherCovariates;
    Ntumor = length(subj_tumor_list);
    Regressors = ones(Ntumor,1);
    Regressors = [SeedSeries,Regressors,OtherCovariatesMatrix];
    Contrast = zeros(1,size(Regressors,2));
    Contrast(1) = 1;
    Predictor = Regressors;
    
    valid_idx = ~isnan(SeedSeries) & all(~isnan(OtherCovariates),2);
    SeedSeries_valid = SeedSeries(valid_idx);
    OtherCovariates_valid = OtherCovariates(valid_idx,:);
    
    OtherCovariatesMatrix = OtherCovariates_valid;
    N_roi = length(SeedSeries_valid);
    Regressors = ones(N_roi,1);
    Regressors = [SeedSeries_valid,Regressors,OtherCovariatesMatrix];
    Contrast = zeros(1,size(Regressors,2));
    Contrast(1) = 1;
    Predictor = Regressors;
    
    for i_GSR = 1:2
        load([work_dir, '/null_distribution/', ...% change here to get the null distribution
            questionnaire_set{i_questionnaire}, GSR_set{i_GSR}, '.mat']);
        t_simul_lh = real(t_simul_lh);
        t_simul_rh = real(t_simul_rh);
        avg_metric_lh = mean(t_simul_lh,1, 'omitnan');
        avg_metric_rh = mean(t_simul_rh,1, 'omitnan');
        sd_metric_lh = std(t_simul_lh, 0, 1, 'omitnan');
        sd_metric_rh = std(t_simul_rh, 0, 1, 'omitnan');
        for i_group = 1:length(group_set)
            if i_GSR == 1
                folder_path = [work_dir, '/mean_FC_mapsS'];
            elseif i_GSR == 2
                folder_path = [work_dir, '/mean_FC_global_mapsS'];
            end
            
            % left hemisphere
            MaskFile = Cfg.MaskFileSurfLH;
            MaskData=gifti(MaskFile);
            MaskData=MaskData.cdata;
            MaskData = double(logical(MaskData));

            % tumor dataset
            all_files = dir(fullfile(folder_path, [group_set{i_group},'_mean_maps/*.gii']));
            is_lh = contains({all_files.name}, '_LH');
            % tumor dataset
            lh_files = fullfile(folder_path, [group_set{i_group},'_mean_maps'], {all_files(is_lh).name}');

            [DependentVolume,~,~,header_lh] = y_ReadAll(lh_files);
            DependentVolume = DependentVolume(:,valid_idx);
            nDimVertex = size(DependentVolume, 1);
            t_surface_lh = zeros(nDimVertex,1);
            CovVariable = []; % no covariate images
            for i_vertex=1:nDimVertex
                if MaskData(i_vertex,1)
                    DependentVariable=DependentVolume(i_vertex,:)';
                    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
                    t_surface_lh(i_vertex,1) = TF_ForContrast;
                end
            end
        
            % right hemisphere
            MaskFile = Cfg.MaskFileSurfRH;
            MaskData=gifti(MaskFile);
            MaskData=MaskData.cdata;
            MaskData = double(logical(MaskData));

            %tumor dataset
            all_files = dir(fullfile(folder_path, [group_set{i_group},'_mean_maps/*.gii']));
            is_rh = contains({all_files.name}, '_RH');
            %tumor dataset
            rh_files = fullfile(folder_path, [group_set{i_group},'_mean_maps'], {all_files(is_rh).name}');

           [DependentVolume,~,~,header_rh] = y_ReadAll(rh_files);
            DependentVolume = DependentVolume(:,valid_idx);
            nDimVertex = size(DependentVolume, 1);
            t_surface_rh = zeros(nDimVertex,1);
            CovVariable = []; % no covariate images
            for i_vertex=1:nDimVertex
                if MaskData(i_vertex,1)
                    DependentVariable=DependentVolume(i_vertex,:)';
                    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(DependentVariable,[Predictor,CovVariable],Contrast,TF_Flag);
                    t_surface_rh(i_vertex,1) = TF_ForContrast;
                end
            end
            
            y_Write(t_surface_lh, header_lh, fullfile(output_dir, ...
                [questionnaire_set{i_questionnaire}, '_', group_set{i_group}, GSR_set{i_GSR}, '_lh.gii']));
            y_Write(t_surface_rh, header_rh, fullfile(output_dir, ...
                [questionnaire_set{i_questionnaire}, '_', group_set{i_group}, GSR_set{i_GSR}, '_rh.gii']));
            
            z_surface_lh = (t_surface_lh - avg_metric_lh')./sd_metric_lh';
            z_surface_rh = (t_surface_rh - avg_metric_rh')./sd_metric_rh';
            
            y_Write(z_surface_lh, header_lh, fullfile(output_dir, ...
                ['z_', questionnaire_set{i_questionnaire}, '_', group_set{i_group}, GSR_set{i_GSR}, '_lh.gii']));
            y_Write(z_surface_rh, header_rh, fullfile(output_dir, ...
                ['z_', questionnaire_set{i_questionnaire}, '_', group_set{i_group}, GSR_set{i_GSR}, '_rh.gii']));
            
        end%group
    end%GSR

 end%questionnaire 