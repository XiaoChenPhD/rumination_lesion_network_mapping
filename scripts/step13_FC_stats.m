% do stats on the FC maps from all the significant rois generated from the
% corrected LNM maps
%
% Xiao Chen 260514
% chenxiaophd@gmail.com

%% initialization

clear; clc;
data_dir = '/xxxxxx/FC_maps_raw';
work_dir = '/xxxxxx/FC_roi_Suzhou';
TemplateDir = '/dpabi_dir/DPABISurf/SurfTemplates';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end

% load participants 
behavior_data = readtable('/xxxxxx/sample_cleaned_merged_v5.csv', 'ReadVariableNames', true);
sub_list = behavior_data.Serial_Number;
Dx = behavior_data.Dx_selected_data;
age = behavior_data.Age;
sex = behavior_data.Sex;
sex = sex -1;
edu = behavior_data.Edu;
HeadMotion = [behavior_data.head_motion_rum, behavior_data.head_motion_dis, ...
              behavior_data.head_motion_rest];
          
scale_data_full = readtable('/xxxxxx/demographic_info_v10.xlsx');
scale_score = scale_data_full(ismember(scale_data_full.Serial_Number, sub_list), :);
scale_name_set = scale_score.Properties.VariableNames(11:end);
vars = varfun(@isnumeric, scale_score, 'OutputFormat', 'uniform');
scale_score{:, vars}(scale_score{:, vars} == 999) = NaN;

HemisphereNameSet = {'lh', 'rh'};

MaskSurfaceRightSet = {'positive_gt3p29_cluster1_n69_mask.func', ...
                       'positive_gt3p29_cluster2_n50_mask.func', ...
                       'positive_gt3p29_cluster3_n33_mask.func'
                            };
                        
GSRSet = {'FunSurfWCF','FunSurfWglobalCF'};
GSRVoluSet = {'FunVoluWCF','FunVoluWglobalCF'};

%% Arrange Files 
% Arrange FC Files
for iGSR = 1:2
    StartingDirName = GSRSet{iGSR};
    StartingDirName_Volume = GSRVoluSet{iGSR};
    for iMask = 1:length(MaskSurfaceRightSet)
        TargetDirRum = [work_dir,'/FC_Maps/MDD_Rum'];
        TargetDirDis = [work_dir,'/FC_Maps/MDD_Dis'];
        if ~exist(TargetDirRum, 'dir'); mkdir(TargetDirRum); end
        if ~exist(TargetDirDis, 'dir'); mkdir(TargetDirDis); end

        mkdir([TargetDirRum,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_lh']); 
        mkdir([TargetDirRum,'/',MaskSurfaceRightSet{iMask}, '_',StartingDirName, '_rh']); 
        mkdir([TargetDirRum,'/', MaskSurfaceRightSet{iMask},'_',StartingDirName, '_subcorticol']); 
        mkdir([TargetDirDis,'/', MaskSurfaceRightSet{iMask},'_',StartingDirName,  '_lh']); 
        mkdir([TargetDirDis,'/', MaskSurfaceRightSet{iMask},'_',StartingDirName, '_rh']); 
        mkdir([TargetDirDis,'/', MaskSurfaceRightSet{iMask}, '_',StartingDirName, '_subcorticol']);

        TargetDirRum = [work_dir,'/FC_Maps/HC_Rum'];
        TargetDirDis = [work_dir,'/FC_Maps/HC_Dis'];
        if ~exist(TargetDirRum, 'dir'); mkdir(TargetDirRum); end
        if ~exist(TargetDirDis, 'dir'); mkdir(TargetDirDis); end
        mkdir([TargetDirRum,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_lh']); 
        mkdir([TargetDirRum,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_rh']); 
        mkdir([TargetDirRum,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_subcorticol']); 
        mkdir([TargetDirDis,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_lh']); 
        mkdir([TargetDirDis,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_rh']); 
        mkdir([TargetDirDis,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_subcorticol']); 
    end
end

  
for iGSR = 1:2
    StartingDirName = GSRSet{iGSR};
    StartingDirName_Volume = GSRVoluSet{iGSR};
    for iMask = 1:length(MaskSurfaceRightSet)
        for i = 1:length(find(Dx == 1))
            % MDD
            TargetDirRum = [work_dir,'/FC_Maps/MDD_Rum'];
            TargetDirDis = [work_dir,'/FC_Maps/MDD_Dis'];
            % rum
            OriginFile = [data_dir,'/ResultsS/FunSurfLH/FC_SeedSurfLHSurfRHVolu_',StartingDirName,'_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirRum,'/', MaskSurfaceRightSet{iMask},'_',StartingDirName, '_lh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

            OriginFile = [data_dir,'/ResultsS/FunSurfRH/FC_SeedSurfLHSurfRHVolu_',StartingDirName, '_', ...
                MaskSurfaceRightSet{iMask},'/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirRum,'/',MaskSurfaceRightSet{iMask}, '_', StartingDirName, '_rh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

        %     OriginFile = []; TargetFile = [];
        %     OriginFile = [data_dir,'/Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/szFC_sub-',sub_list{i},'.nii'];
        %     TargetFile = [TargetDirRum,'/subcorticol/szFC_sub-',sub_list{i},'.nii'];
        %     copyfile(OriginFile,TargetFile);


            OriginFile = [data_dir,'/S2_ResultsS/FunSurfLH/FC_SeedSurfLHSurfRHVolu_',StartingDirName, '_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirDis,'/', MaskSurfaceRightSet{iMask},'_', StartingDirName, '_lh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

            OriginFile = [data_dir,'/S2_ResultsS/FunSurfRH/FC_SeedSurfLHSurfRHVolu_',StartingDirName,'_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirDis,'/', MaskSurfaceRightSet{iMask},'_', StartingDirName, '_rh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

        %     OriginFile = []; TargetFile = [];
        %     OriginFile = [data_dir,'/Task/S2_Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
        %     TargetFile = [TargetDirDis,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
        %     copyfile(OriginFile,TargetFile);
        end

        % HC
        TargetDirRum = [work_dir,'/FC_Maps/HC_Rum'];
        TargetDirDis = [work_dir,'/FC_Maps/HC_Dis'];
        % rum
        for i = length(find(Dx == 1))+1:length(sub_list)

            OriginFile = [data_dir,'/ResultsS/FunSurfLH/FC_SeedSurfLHSurfRHVolu_',StartingDirName, '_', ...
                MaskSurfaceRightSet{iMask},'/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirRum,'/', MaskSurfaceRightSet{iMask},'_', StartingDirName, '_lh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

            OriginFile = [data_dir,'/ResultsS/FunSurfRH/FC_SeedSurfLHSurfRHVolu_',StartingDirName, '_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirRum,'/', MaskSurfaceRightSet{iMask}, '_',StartingDirName, '_rh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        %     OriginFile = []; TargetFile = [];
        %     OriginFile = [data_dir,'/Task/Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
        %     TargetFile = [TargetDirRum,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
        %     copyfile(OriginFile,TargetFile);

            OriginFile = [data_dir,'/S2_ResultsS/FunSurfLH/FC_SeedSurfLHSurfRHVolu_',StartingDirName,'_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirDis,'/',MaskSurfaceRightSet{iMask},'_', StartingDirName, '_lh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

            OriginFile = [data_dir,'/S2_ResultsS/FunSurfRH/FC_SeedSurfLHSurfRHVolu_',StartingDirName, '_', ...
                MaskSurfaceRightSet{iMask}, '/szFC_sub-',sub_list{i},'.func.gii'];
            TargetFile = [TargetDirDis,'/',MaskSurfaceRightSet{iMask},'_',StartingDirName, '_rh/szFC_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);

        %     OriginFile = []; TargetFile = [];
        %     OriginFile = [data_dir,'/S2_ResultsS/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/szFC_sub-',sub_list{i},'.nii'];
        %     TargetFile = [TargetDirDis,'/subcorticol/szFC_sub-',sub_list{i},'.nii'];
        %     copyfile(OriginFile,TargetFile);
        end
    end
end

%% do paired t test on surface

PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

for iGSR = 1:2
    StartingDirName = GSRSet{iGSR};
    StartingDirName_Volume = GSRVoluSet{iGSR};
    for iMask = 1:length(MaskSurfaceRightSet)
        TopOutput_Dir = [work_dir,'/stats/',MaskSurfaceRightSet{iMask},'_', StartingDirName];
        if ~exist(TopOutput_Dir, 'dir'); mkdir(TopOutput_Dir); end

        % MDD
        OtherCovariates = {};
        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);

        for iHem = 1:2
            DependentDirs{1,1} = [work_dir,'/FC_Maps/MDD_Rum/',MaskSurfaceRightSet{iMask},'_', ...
                StartingDirName,'_', HemisphereNameSet{iHem}];
            DependentDirs{2,1} = [work_dir,'/FC_Maps/MDD_Dis/',MaskSurfaceRightSet{iMask},'_', ...
                StartingDirName, '_', HemisphereNameSet{iHem}];

            PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white.surf.gii'];
            PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white_avg.area.gii'];
            MaskFile =  [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_cortex.label.gii'];

            OutputName = [TopOutput_Dir,'/PT_MDD_',HemisphereNameSet{iHem}];
            [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
        end


        % HC
        OtherCovariates = {};
        OtherCovariates{1,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        for iHem = 1:2
            DependentDirs{1,1} = [work_dir,'/FC_Maps/HC_Rum/',MaskSurfaceRightSet{iMask},'_', ...
                StartingDirName,'_', HemisphereNameSet{iHem}];
            DependentDirs{2,1} = [work_dir,'/FC_Maps/HC_Dis/',MaskSurfaceRightSet{iMask},'_', ...
                StartingDirName, '_', HemisphereNameSet{iHem}];
            PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white.surf.gii'];
            PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white_avg.area.gii'];
            MaskFile =  [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_cortex.label.gii']; 
            OutputName = [TopOutput_Dir,'/PT_HC_',HemisphereNameSet{iHem}];
            [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
        end
    end
end


%%  do 7 network analysis with Schaefer 400 parcelations
TemplateDir = '/xxxxxx/atlas';
Schaefer_400_17_lh = [TemplateDir, '/fsaverage5_lh_Schaefer2018_400Parcels_17Networks_order.label.gii'];
Schaefer_400_17_rh = [TemplateDir, '/fsaverage5_rh_Schaefer2018_400Parcels_17Networks_order.label.gii'];
dpabiSurfTemplateDir = '/dpabi_dir/DPABISurf/SurfTemplates';
HemisphereNameSet = {'lh', 'rh'};
ConditionSet = {'MDD_Rum', 'MDD_Dis', 'HC_Rum', 'HC_Dis'};

%% extract the FC matrix
mask_name_set = {'cluster1', 'cluster2', 'cluster3'};
OutputDir = [work_dir, '/FC_network_extracted'];
if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
for iGSR = 1:2
    StartingDirName = GSRSet{iGSR};
    StartingDirName_Volume = GSRVoluSet{iGSR};
    for iCondition = 1:length(ConditionSet)
        mkdir([OutputDir,filesep, ConditionSet{iCondition}]);
        for iMask = 1:length(MaskSurfaceRightSet)
            AllVolume = y_ReadAll([work_dir, '/FC_Maps/', ConditionSet{iCondition}, '/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh']);
            [ROISignalsSurfLH] = y_ExtractROISignal_Surf(AllVolume, {Schaefer_400_17_lh}, ...
                                        [OutputDir,filesep,ConditionSet{iCondition}, filesep, mask_name_set{iMask},'_', StartingDirName,'_lh'], ...
                                        '', 1);

            AllVolume = y_ReadAll([work_dir, '/FC_Maps/', ConditionSet{iCondition}, '/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh']);                      
            [ROISignalsSurfRH] = y_ExtractROISignal_Surf(AllVolume, {Schaefer_400_17_rh}, ...
                                        [OutputDir,filesep,ConditionSet{iCondition}, filesep, mask_name_set{iMask},'_', StartingDirName,'_rh'], ...
                                        '', 1);
        end
    end
end

%% get Network FC
load('~/rumination_lesion_network_mapping/data/DPABISurf_Schaefer2018_400_Tian2020_54_Info.mat');

% network_FC: 1 x 4 cells, 4 conditions in ConditionSet,
% each cell has n x 10 matrixs, 10 networks defined as the abovementrioned
% .mat file
% no_GSR
mask_name_set = {'cluster1', 'cluster2', 'cluster3'};
TopOutputDir = [work_dir, '/FC_network_extracted'];

% no GSR
network_FC_all = struct();
for iMask = 1:length(MaskSurfaceRightSet)
    network_FC = {};
    for iCondition = 1:length(ConditionSet)
        load([TopOutputDir, '/', ConditionSet{iCondition}, '/ROISignals_',mask_name_set{iMask}, '_FunSurfWCF_lh.mat']);
        ROISignals_lh = ROISignals;
        load([TopOutputDir, '/', ConditionSet{iCondition}, '/ROISignals_',mask_name_set{iMask}, '_FunSurfWCF_rh.mat']);
        ROISignals_rh = ROISignals;
%         load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_subcortical.mat']);
%         ROISignals_subcortical = ROISignals;
        ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh];
%         ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh,ROISignals_subcortical];
        ROISignals_full_384{iCondition} = ROISignals_full{iCondition}(:, ROIIndex_181_Schafer_Tian(1,1:384));
        network_label = DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(1:384,1)';
        for iNetwork = unique(network_label)
            network_FC{iCondition}(:,iNetwork) = mean(ROISignals_full_384{iCondition}(:, find(network_label == iNetwork)), 2);
        end
    end
    network_FC_all.(mask_name_set{iMask}) = network_FC;
end

% GSR
network_FC_all_GSR = struct();
for iMask = 1:length(MaskSurfaceRightSet)
    network_FC = {};
    for iCondition = 1:length(ConditionSet)
        load([TopOutputDir, '/', ConditionSet{iCondition}, '/ROISignals_',mask_name_set{iMask}, '_FunSurfWglobalCF_lh.mat']);
        ROISignals_lh = ROISignals;
        load([TopOutputDir, '/', ConditionSet{iCondition}, '/ROISignals_',mask_name_set{iMask}, '_FunSurfWglobalCF_rh.mat']);
        ROISignals_rh = ROISignals;
%         load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_subcortical.mat']);
%         ROISignals_subcortical = ROISignals;
        ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh];
%         ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh,ROISignals_subcortical];
        ROISignals_full_384{iCondition} = ROISignals_full{iCondition}(:, ROIIndex_181_Schafer_Tian(1,1:384));
        network_label = DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork(1:384,1)';
        for iNetwork = unique(network_label)
            network_FC{iCondition}(:,iNetwork) = mean(ROISignals_full_384{iCondition}(:, find(network_label == iNetwork)), 2);
        end
    end
    network_FC_all_GSR.(mask_name_set{iMask}) = network_FC;
end

%% do mixed effect analysis
mask_name_set = {'cluster1', 'cluster2', 'cluster3'};
TopOutputDir = [work_dir, '/stats_FC_network'];
if ~exist(TopOutputDir, 'dir'); mkdir(TopOutputDir); end

for iMask = 1:length(MaskSurfaceRightSet)
    for iGSR = 1:2
        if iGSR == 1
            network_FC = network_FC_all.(mask_name_set{iMask});
        else
            network_FC = network_FC_all_GSR.(mask_name_set{iMask});
        end

        MDD_Rum_n = size(network_FC{1},1);
        MDD_Dis_n = size(network_FC{2},1);
        HC_Rum_n = size(network_FC{3},1);
        HC_Dis_n = size(network_FC{4},1);

        BetweenSubjectFactor1 = ones(MDD_Rum_n+MDD_Dis_n,1);
        BetweenSubjectFactor2 = -1*ones(HC_Rum_n+HC_Dis_n,1);
        BetweenSubjectFactor = [BetweenSubjectFactor1;BetweenSubjectFactor2];
        WithinSubjectFactor = [ones(MDD_Rum_n,1);-1*ones(MDD_Dis_n,1);ones(HC_Rum_n,1);-1*ones(HC_Dis_n,1)];
        nSubject = [MDD_Rum_n MDD_Dis_n HC_Rum_n HC_Dis_n];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        SubjectRegressorsAll = [];
        for i=1:4
            SubjectRegressorsAll=[SubjectRegressorsAll;(ceil(i/2)-1)*10000+([1:nSubject(i)]')];
        end
        SubIndex=unique(SubjectRegressorsAll);
        nSub=length(SubIndex);
        SubjectRegressors=[];
        for i=1:nSub
            SubjectRegressors(:,i) = zeros(size(SubjectRegressorsAll));
            SubjectRegressors(SubjectRegressorsAll==SubIndex(i),i) = 1;
        end
        Interaction = WithinSubjectFactor.*BetweenSubjectFactor;
        HeadMotion_Regressor = [OtherCovariates{1};OtherCovariates{2};OtherCovariates{3};OtherCovariates{4}];
        AllCov = [WithinSubjectFactor,Interaction,SubjectRegressors,HeadMotion_Regressor];
        Contrast = zeros(1,size(AllCov,2));
        TF_Flag = 'T';

        % condition effect
        Contrast(1)=1;
        for iNetwork = unique(network_label)
            y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
            Stats_Mixed_Condition{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
            Stats_Mixed_Condition{iNetwork,2} = TF_ForContrast;
            Stats_Mixed_Condition{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
            Stats_Mixed_Condition{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
            Stats_Mixed_Condition{iNetwork,5} = Cohen_f2;
            %p_raw(iNetwork) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
        end

        % interaction effect
        Contrast(1)=0;
        Contrast(2)=1;
        for iNetwork = unique(network_label)
            y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
            [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
            Stats_Mixed_Interaction{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
            Stats_Mixed_Interaction{iNetwork,2} = TF_ForContrast;
            Stats_Mixed_Interaction{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
            Stats_Mixed_Interaction{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
            Stats_Mixed_Interaction{iNetwork,5} = Cohen_f2;
            p_raw(iNetwork) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
        end
        % % fdr correction
        % p_raw = p_raw(1:9);
        % p_fdr = mafdr(p_raw, 'BHFDR',true);
        % FDRMsk=y_FDR_Vector(p_raw', 0.05);
        
        if iGSR == 1
            save([TopOutputDir,'/', mask_name_set{iMask}, '.mat'], 'Stats_Mixed_Interaction', 'Stats_Mixed_Condition');
        else
            save([TopOutputDir,'/', mask_name_set{iMask}, '_GSR.mat'], 'Stats_Mixed_Interaction', 'Stats_Mixed_Condition');
        end
        
    end
end

%% scale corr

mask_name_set = {'cluster1', 'cluster2', 'cluster3'};
TopOutputDir = [work_dir, '/stats_FC_network'];
if ~exist(TopOutputDir, 'dir'); mkdir(TopOutputDir); end

scale_data_full = readtable('/xxxxxxx/demographic_info_v10.xlsx');
scale_score = scale_data_full(ismember(scale_data_full.Serial_Number, sub_list), :);
scale_name_set = scale_score.Properties.VariableNames(11:end);
vars = varfun(@isnumeric, scale_score, 'OutputFormat', 'uniform');
scale_score{:, vars}(scale_score{:, vars} == 999) = NaN;

% construct MRI Phenotype
% for each network FCs (cluster GSR vs. no GSR), construct a cell to store all the values
% sig_cluster1, sig_cluster2, sigcluster3, GSR, no GSR
% a 1 x 8 cell, each cell contains a n x 9 matrix:
% cell 1-4: MDD_Rum, MDD_Dis, HC_Rum, HC_Dis
% cell 5-6: MDD_Rum-MDD_Dis, HC_Rum-HC_Dis
% cell 7-8: (MDD_Rum-MDD_Dis)/MDD_Dis, (HC_Rum-HC_Dis)/HC_Dis

MRI_label = YeoSCNetwork_Label(1:9,2); %j
condition_label = {'MDD_Rum', 'MDD_Dis', 'HC_Rum', 'HC_Dis', ...
                   'MDD_delta', 'HC_delta', ...
                   'MDD_relative_delta', 'HC_relative_delta'}; %i

for iMask = 1:length(MaskSurfaceRightSet)
    for iGSR = 1:2
        if iGSR == 1
            MRI_Phenotype = network_FC_all.(mask_name_set{iMask});
        else
            MRI_Phenotype = network_FC_all_GSR.(mask_name_set{iMask});
        end

        MRI_Phenotype{5} = MRI_Phenotype{1} - MRI_Phenotype{2};
        MRI_Phenotype{6} = MRI_Phenotype{3} - MRI_Phenotype{4};
        MRI_Phenotype{7} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
        MRI_Phenotype{8} = (MRI_Phenotype{3} - MRI_Phenotype{4})./MRI_Phenotype{4};

        % do scale corr analyses, explorative
        % MDD patients
        stats_corr = [];
        temp_MDD = [];
        for i = [1,2,5,7] %select MDD
            count = 0;
            for j = 1:size(MRI_Phenotype{1},2)
                for k = 11:size(scale_score,2) % from columne 11 to end
                    y = MRI_Phenotype{i}(:,j);
                    x = scale_score{Dx == 1, k};
                    [r, p] = corr(x,y, 'Rows', 'complete');
                    if p < 0.005
                        count = count + 1;
                        eval(['stats_corr.',condition_label{i},'{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k-10}]']);
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

        % HC
        temp_HC = [];
        for i = [3,4,6,8] % select HC
            count = 0;
            for j = 1:size(MRI_Phenotype{1},2)
                for k = 11:size(scale_score,2) % from columne 11 to end

                    y = MRI_Phenotype{i}(:,j);
                    x = scale_score{Dx == 2, k};
                    [r, p] = corr(x,y, 'Rows', 'complete');
                    if p < 0.005
                        count = count + 1;
                        eval(['stats_corr.',condition_label{i},'{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k-10}]']);
                        eval(['stats_corr.',condition_label{i},'{count, 2} = r']);
                        eval(['stats_corr.',condition_label{i},'{count, 3} = p']);
                        %store the significant condition ids
                        temp_HC(end+1,1) = i;
                        temp_HC(end,2) = j;
                        temp_HC(end,3) = k;
                     end
                end
            end
        end
        
        if iGSR == 1
            save([TopOutputDir,'/',mask_name_set{iMask}, '_corr.mat'], 'stats_corr');
        else
            save([TopOutputDir,'/',mask_name_set{iMask}, '_GSR_corr.mat'], 'stats_corr');
        end
    end
end


%% legacy code： do mixed effect analysis on surface
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

for iGSR = 1:2
    StartingDirName = GSRSet{iGSR};
    StartingDirName_Volume = GSRVoluSet{iGSR};
    for iMask = 1:length(MaskSurfaceRightSet)

        % left hemisphere
        PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_lh_white.surf.gii'];
        PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_lh_white_avg.area.gii'];
        MaskFile =  [TemplateDir,'/fsaverage5_lh_cortex.label.gii'];

        DependentDir{1,1} =  [work_dir,'/FC_Maps/MDD_Rum/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh'];
        DependentDir{2,1} = [work_dir,'/FC_Maps/MDD_Dis/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh'];
        DependentDir{3,1} = [work_dir,'/FC_Maps/HC_Rum/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh'];
        DependentDir{4,1} = [work_dir,'/FC_Maps/HC_Dis/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh'];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        OutputDir = [work_dir,'/stats/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_lh'];
        if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
        OutputName = [OutputDir,'/FC_stats'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);

        % right hemisphere
        PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_rh_white.surf.gii'];
        PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_rh_white_avg.area.gii'];
        MaskFile =  [TemplateDir,'/fsaverage5_rh_cortex.label.gii'];

        DependentDir{1,1} =  [work_dir,'/FC_Maps/MDD_Rum/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh'];
        DependentDir{2,1} = [work_dir,'/FC_Maps/MDD_Dis/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh'];
        DependentDir{3,1} = [work_dir,'/FC_Maps/HC_Rum/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh'];
        DependentDir{4,1} = [work_dir,'/FC_Maps/HC_Dis/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh'];

        OutputDir = [work_dir,'/stats/',MaskSurfaceRightSet{iMask},'_', StartingDirName,'_rh'];
        if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
        OutputName = [OutputDir,'/FC_stats'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);

        % % volume
        % PALMSettings = [];
        % PALMSettings.nPerm = 5000;
        % PALMSettings.ClusterInference=0;
        % PALMSettings.ClusterFormingThreshold=2.3;
        % PALMSettings.TFCE=1;
        % PALMSettings.FDR=0;
        % PALMSettings.TwoTailed=1;
        % PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
        % PALMSettings.SavePermutations = 0;
        % 
        % MaskFile =  [data_dir,'/Task/Masks/AllResampled_BrainMask_05_91x109x91.nii'];
        % 
        % DependentDir{1,1} =  [work_dir,'/MDD_Rum/subcorticol'];
        % DependentDir{2,1} = [work_dir,'/MDD_Dis/subcorticol'];
        % DependentDir{3,1} = [work_dir,'/HC_Rum/subcorticol'];
        % DependentDir{4,1} = [work_dir,'/HC_Dis/subcorticol'];
        % 
        % OutputDir = [work_dir,'/FC_analysis/stats/subcorticol'];
        % if ~exist(OutputDir); mkdir(OutputDir); end
        % OutputName = [OutputDir,'/FC_stats'];
        % y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
    end
end

