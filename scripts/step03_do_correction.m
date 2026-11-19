% Using permutation test to do the correction.
% according to Cash et al., 2023. Nature Mental Health
%
% generate random rois, and calculate the Fisher's r-to-z
% transformed FC
% 
% accelarating the calculation
% 251107
%
% v2: further accelarating calculating using matrix 
% 251108
%
% Wrapped up by Xiao Chen 261119
%
% Xiao Chen
% chenxiaophd@gmail.com

%% initialization
clear;clc;

addpath('~/rumination_lesion_network_mapping/scripts');
work_dir = 'xxxxxxxxxxxxx';
DIRECT_mask = '~/rumination_lesion_network_mapping/data/AllResampled_BrainMask_05_91x109x91.nii';
% the DIRECT subj list
SubjectListFile = 'xxxxxx/sub_list.txt';

load('xxxxxxx/Cfg_3_Processing.mat');
Cfg.WorkingDir = 'xxxxxxxx';
Cfg.MaskFileVolu = '~/rumination_lesion_network_mapping/data/AllResampled_BrainMask_05_91x109x91.nii';
Cfg.SubjectID = importdata(SubjectListFile);

GSRSet = {'FunSurfWCF','FunSurfWglobalCF'};
GSRVoluSet = {'FunVoluWCF','FunVoluWglobalCF'};

[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));
ProgramPath=fullfile(DPABIPath, 'DPARSF');
Cfg.MaskFileSurfLH = [DPABIPath, '/DPABISurf/SurfTemplates/fsaverage5_lh_cortex.label.gii'];
Cfg.MaskFileSurfRH = [DPABIPath, '/DPABISurf/SurfTemplates/fsaverage5_rh_cortex.label.gii'];

% read in the lesion subjects' list
demog_info = readtable('/xxxxxx/scale_data_raw.xlsx', 'sheet', 1);
subj_tumor_list = demog_info{:,1};
subj_tumor_list(1,:) =[]; %remove the example row
subj_tumor_list = cellfun(@lower, subj_tumor_list, 'UniformOutput', false); %lower_case

%% open the parallel sessions

poolobj = gcp('nocreate');
if isempty(poolobj) || poolobj.NumWorkers ~= 30
    delete(poolobj);
    parpool(30);
end

%% ramdomly pick rois
% get the necessary info
TotalPerm = 1000;       % total number of permutations
BatchSize = 20;         % number of permutations per batch
NumBatches = TotalPerm / BatchSize;

output_dir = fullfile(work_dir, 'xxxxxx');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

Ntumor = length(subj_tumor_list);
Radius = 5;
RefFile = DIRECT_mask;

[mask_data_surf, ~] = y_ReadAll(Cfg.MaskFileSurfLH);
vertex_length = length(mask_data_surf);

V = spm_vol(DIRECT_mask);
mask_data = spm_read_vols(V);
idx = find(mask_data == 1);

% dir files
Nsub = length(Cfg.SubjectID);
VolFile  = cell(Nsub, 2); % (subject, iGSR)
SurfLFile = cell(Nsub, 2);
SurfRFile = cell(Nsub, 2);

for iGSR = 1:2
    Cfg.StartingDirName        = GSRSet{iGSR};
    Cfg.StartingDirName_Volume = GSRVoluSet{iGSR};

    for i = 1:Nsub
        subj = Cfg.SubjectID{i};

        d = dir(fullfile(Cfg.WorkingDir, Cfg.StartingDirName_Volume, ...
                         subj, 'Filtered_4DVolume.nii'));
        VolFile{i, iGSR} = fullfile(d(1).folder, d(1).name);

        dL = dir(fullfile(Cfg.WorkingDir, Cfg.StartingDirName, subj, ...
                          '*fsaverage5_hemi-L*.func.gii'));
        SurfLFile{i, iGSR} = fullfile(dL(1).folder, dL(1).name);

        dR = dir(fullfile(Cfg.WorkingDir, Cfg.StartingDirName, subj, ...
                          '*fsaverage5_hemi-R*.func.gii'));
        SurfRFile{i, iGSR} = fullfile(dR(1).folder, dR(1).name);
    end
end

for ibatch = 1:NumBatches
    fprintf('Running batch %d / %d (%d–%d)...\n', ...
        ibatch, NumBatches, (ibatch-1)*BatchSize+1, ibatch*BatchSize);
    
    % generate all the random spheres, time consuming...
    SphereIdx = cell(BatchSize, Ntumor); 
    ReDrawCount = zeros(BatchSize, Ntumor);
    max_attempts = 50;
    for ip = 1:BatchSize
        for it = 1:Ntumor
            validSphere = false;
            attempt = 0;
            while ~validSphere && attempt < max_attempts
                attempt = attempt + 1;
                sel = idx(randi(numel(idx)));
                [x, y, z] = ind2sub(size(mask_data), sel);
                MNIcoord = V.mat * [x; y; z; 1];
                Center = MNIcoord(1:3)';
                [SphereData, SHeader] = y_Sphere(Center, Radius, RefFile);
                SphereData = SphereData .* mask_data;

                % check if the sphere is located 
                % in the mask
                voxIdx = find(SphereData);
                if isempty(voxIdx)
                    continue;
                end
                ratio = numel(voxIdx)*8 / ((4/3)*pi*(Radius^3));
                if ratio < 0.3
                    continue;
                end
                %y_Write(SphereData,SHeader,'test.nii')
                validSphere = true;
            end
            SphereIdx{ip, it} = find(SphereData); 
            if ~validSphere
                warning('failed to get valid spheres after %d attempts, using the last location', ...
            max_attempts);
            end
            ReDrawCount(ip, it) = attempt - 1;
        end
    end
  
    PermMean_lh     = zeros(BatchSize, Ntumor, vertex_length, 'single');
    PermMean_rh     = zeros(BatchSize, Ntumor, vertex_length, 'single');
    PermMean_lh_GSR = zeros(BatchSize, Ntumor, vertex_length, 'single');
    PermMean_rh_GSR = zeros(BatchSize, Ntumor, vertex_length, 'single');

    % main loop
    for iGSR = 1:2
        fprintf('Processing GSR condition %d...\n', iGSR);
        subj_lh  = cell(Nsub, 1);%each cell: BatchSize × Ntumor × vertex
        subj_rh  = cell(Nsub, 1);
        PermMean_lh_cur = zeros(BatchSize, Ntumor, vertex_length, 'single');
        PermMean_rh_cur = zeros(BatchSize, Ntumor, vertex_length, 'single');

        parfor iSub = 1:Nsub
            [Data4D, ~,~,Header] = y_ReadAll(VolFile{iSub, iGSR});
            T = size(Data4D, 4);
            Data2D = reshape(Data4D, [], T)';  % T × Nvoxel
            fileL = SurfLFile{iSub, iGSR};
            fileR = SurfRFile{iSub, iGSR};

            % 1. seed matrix
            Nseed = BatchSize * Ntumor;
            SeedMat = zeros(T, Nseed, 'single');
            k = 0;
            for it = 1:Ntumor
                for ip = 1:BatchSize
                    k = k + 1;
                    voxIdx = SphereIdx{ip, it};
                    SeedMat(:, k) = mean(Data2D(:, voxIdx), 2, 'omitnan');
                end
            end

            % 2. read in the surface data
            [BOLD_lh, ~, ~, GHeader_lh] = y_ReadAll(fileL);   % nVert × T
            [BOLD_rh, ~, ~, GHeader_rh] = y_ReadAll(fileR);

            % 3. LH / RH 
            z_lh_sub = fast_sca_surf_batch(BOLD_lh, Cfg.MaskFileSurfLH, SeedMat);
            z_rh_sub = fast_sca_surf_batch(BOLD_rh, Cfg.MaskFileSurfRH, SeedMat);

            % 4. reshape
            z_lh_sub = reshape(z_lh_sub, BatchSize, Ntumor, vertex_length);
            z_rh_sub = reshape(z_rh_sub, BatchSize, Ntumor, vertex_length);

            % 5. reduction
            PermMean_lh_cur = PermMean_lh_cur + z_lh_sub;
            PermMean_rh_cur = PermMean_rh_cur + z_rh_sub;
        end
    
        PermMean_lh_cur = PermMean_lh_cur / Nsub;
        PermMean_rh_cur = PermMean_rh_cur / Nsub;

        if iGSR == 1
            PermMean_lh = PermMean_lh_cur;
            PermMean_rh = PermMean_rh_cur;
        else
            PermMean_lh_GSR = PermMean_lh_cur;
            PermMean_rh_GSR = PermMean_rh_cur;
        end
    end
    save_filename = sprintf('PermMean_batch%03d.mat', ibatch);
    save(fullfile(output_dir, save_filename), ...
        'PermMean_lh', 'PermMean_rh', ...
        'PermMean_lh_GSR', 'PermMean_rh_GSR', ...
        'SphereIdx', 'ReDrawCount', '-v7.3');

    fprintf('Saved %s\n', fullfile(output_dir, save_filename));

    % (optional) clear large variables to free memory
    clear PermMean_lh PermMean_rh PermMean_lh_GSR PermMean_rh_GSR SphereIdx ReDrawCount
end

