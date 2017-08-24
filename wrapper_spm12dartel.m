%% SPM12 Dartel using just MPRAGE or both MPRAGE & MBW
% Created by Kevin Tan on Jun 29, 2017 (some code adopted from Bob Spunt)

% Instructions, agorithmic description & edit history (please read!!):
%   https://github.com/scanUCLA/spm12-dartel

% Last revision: 2 Aug 2017 - Kevin Tan
%% User-editable Parameters

% Path/directory/name information
owd = '/u/project/sanscn/kevmtan/scripts/rawTestData/PTSD';  % base study directory
codeDir = '/u/project/sanscn/kevmtan/scripts/SPM12_DARTEL'; % where code lives
batchDir = '/u/project/sanscn/kevmtan/scripts/SPM12_DARTEL/PTSD_170802'; % dir in which to save batch scripts & predartel subject status + workspace

% pattern for finding subject folders (use wildcards)
subID = 'VET*';

% Subjects to do/skip, example: {'sub001' 'sub002'}
subNam = {}; % do which subjects? (leave empty to do all)
skipSub = {}; % skip which subjects? (leave empty to do all)

% Funcitonal image info
runID = 'BOLD_*'; % pattern for finding functional run folders (use wildcards)
funcID ='BOLD_'; % first character(s) in your functional NIFTI files? (do NOT use wildcards)

% Are you using 4D NIFTI functional files?
fourDnii = true; % true=all run volumes in one NIFTI files, false=run volumes in separate NIFTI files

% Structural image info (pattern for finding struct folders, assumes first characters of struct folder and struct filename are same)
structID{1} = 'SAG_MPRAGE*'; % MPRAGE if you have it, MBW if no MPRAGE (use wildcards)
structID{2} = 'Matched_Bandwidth_HiRes*'; % MBW if you have it in addition to MPRAGE (use wildcards)

% Two structural scans or just one?
twoStructs = true; % true = 2 structural scans, false = 1 strutural scan

% Path of TPM tissues in your SPM directory
tpmPath = '/u/project/CCN/apps/spm12/tpm';

% Voxel size for resampling (use AFNI's dicom_hdr on the functional & structural DICOM files and use the "slice thickness")
fVoxSize = [3 3 3]; % functionals (non-multiband usually [3 3 3], multiband usually [2 2 2])
sVoxSize = [1 1 1]; % for structID{1} (MPRAGE usually [1 1 1], MBW usually [3 3 3])

% smoothing kernel for functionals  (mm isotropic)
FWHM = 8;

% Number of workers (threads) matlab should use
nWorkers = maxNumCompThreads; % specify integer if desired

% EXECUTE (1) or just make matlabbatches (0)
execPreDartel = 0;
execDartel = 0;

%% Setup subjects

% Make batch folder
try
    mkdir(batchDir);
catch
end

diary([batchDir '/preDARTEL_log_' datestr(now,'yyyymmdd_HHMM') '.txt']);

% Find subject directories
if isempty(subNam)
    d = dir([owd '/' subID]);
    for ii = 1:length(d)
        subNam{ii} = d(ii).name;
        fprintf('Adding %s\n', subNam{ii})
    end
end
numSubs = length(subNam);
cd(codeDir);

% Prepare status struct
runStatus = struct('subNam',[],'status',[],'error',[]);
runStatus(numSubs).subNam = [];
runStatus(numSubs).status = [];
runStatus(numSubs).error = [];

%% Run Pre-dartel (explicitly parallelized per subject)

% Determine number of parallel workers
nWorkers = min(numSubs, nWorkers);
parpool('local', nWorkers);

% Parfor loop to run explicity parallelized pre-dartel across subs
parfor i = 1:numSubs
    % Pre-allocate subject in runStatus struct
    runStatus(i).subNam = subNam{i};
    
    % Cross-check subject with run/skip list
    if ismember(subNam{i}, skipSub)
        runStatus(i).status = 0;
        runStatus(i).error = 'Subject in exclusion list';
        disp(['Skipping subject ' subNam{i} ', is in exclusion list']);
        continue
    else % Run subject
        disp(['Running subject ' subNam{i}]);
        [runStatus(i).status, runStatus(i).error, allfuncs{i}, allt1{i}, allmt1{i},...
            allrc1{i}, allrc2{i}, allu_rc1{i}, allc1{i}, allc2{i}, allc3{i}] =...
            run_spm12dartel(subNam{i}, owd, codeDir, batchDir, runID, funcID,...
            structID, twoStructs, fourDnii, execPreDartel);
        if runStatus(i).status == 1
            disp(['subject ' subNam{i} ' successful']);
        else
            runStatus(i).status = 0;
            disp([runStatus(i).error ' for ' subNam{i}]);
        end
    end
end
delete(gcp('nocreate'));

% Save stuff
date = datestr(now,'yyyymmdd_HHMM');
filename = [batchDir '/runStatus_' date '.mat'];
save(filename,'runStatus');
filename = [batchDir '/forDartel_workspace_' date '.mat']; % Use this to re-do "run dartel" if it fails
save(filename);
diary off

%% Create DARTEL matlabbatch

diary([batchDir '/DARTEL_log_' datestr(now,'yyyymmdd_HHMM') '.txt']);

% Load SPM
spm('defaults','fmri'); spm_jobman('initcfg');

% Set max threads for implicit multithreading
%lastNumThreads = maxNumCompThreads(nWorkers*2);

% Subjects with no errors in pre-dartel
inds = find([runStatus.status]);

% Create Segmentation matlabbatch (implicitly parallelized in SPM)
for s = 1:length(inds)
    matlabbatch{1}.spm.spatial.preproc.channel.vols{s,1} = allt1{inds(s)};
end
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = cellstr([tpmPath '/TPM.nii,1']);
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = cellstr([tpmPath '/TPM.nii,2']);
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = cellstr([tpmPath '/TPM.nii,3']);
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = cellstr([tpmPath '/TPM.nii,4']);
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = cellstr([tpmPath '/TPM.nii,5']);
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = cellstr([tpmPath '/TPM.nii,6']);
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

% DARTEL: create templates
for s = 1:length(inds)
    matlabbatch{2}.spm.tools.dartel.warp.images{1,1}{s,1} = allrc1{inds(s)};
    matlabbatch{2}.spm.tools.dartel.warp.images{1,2}{s,1} = allrc2{inds(s)};
end
matlabbatch{2}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{2}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.its = 3;

% DARTEL: normalize & smooth functional images to MNI
matlabbatch{3}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{3}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{3}.spm.tools.dartel.mni_norm.data.subj(s).images = allfuncs{inds(s)};
end                                             
matlabbatch{3}.spm.tools.dartel.mni_norm.vox = fVoxSize;
matlabbatch{3}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{3}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{3}.spm.tools.dartel.mni_norm.fwhm = [FWHM FWHM FWHM];

% DARTEL: normalize bias-corrected MPRAGE to MNI
matlabbatch{4}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{4}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{4}.spm.tools.dartel.mni_norm.data.subj(s).images{1} = allmt1{inds(s)};
end                                            
matlabbatch{4}.spm.tools.dartel.mni_norm.vox = sVoxSize;
matlabbatch{4}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{4}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{4}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

% DARTEL: normalize grey matter segmentation to MNI (for Conn toolbox or other use of segmentations)
matlabbatch{5}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{5}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{5}.spm.tools.dartel.mni_norm.data.subj(s).images{1} = allc1{inds(s)};
end                                            
matlabbatch{5}.spm.tools.dartel.mni_norm.vox = sVoxSize;
matlabbatch{5}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{5}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{5}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

% DARTEL: normalize white matter segmentation to MNI (for Conn toolbox or other use of segmentations)
matlabbatch{6}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{6}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{6}.spm.tools.dartel.mni_norm.data.subj(s).images{1} = allc2{inds(s)};
end                                            
matlabbatch{6}.spm.tools.dartel.mni_norm.vox = sVoxSize;
matlabbatch{6}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{6}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{6}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

% DARTEL: normalize CSF segmentation to MNI (for Conn toolbox or other use of segmentations)
matlabbatch{7}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{7}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{7}.spm.tools.dartel.mni_norm.data.subj(s).images{1} = allc3{inds(s)};
end                                            
matlabbatch{7}.spm.tools.dartel.mni_norm.vox = sVoxSize;
matlabbatch{7}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{7}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{7}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

% DARTEL: normalize functional images to MNI no smoothing
matlabbatch{8}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
for s = 1:length(inds)
    matlabbatch{8}.spm.tools.dartel.mni_norm.data.subj(s).flowfield{1} = allu_rc1{inds(s)};
    matlabbatch{8}.spm.tools.dartel.mni_norm.data.subj(s).images = allfuncs{inds(s)};
end                                             
matlabbatch{8}.spm.tools.dartel.mni_norm.vox = fVoxSize;
matlabbatch{8}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{8}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{8}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

%% Save & run DARTEL matlabbatch

% save matlabbatch struct in output folder
time_stamp = datestr(now,'yyyymmdd_HHMM');
filename = [batchDir '/DARTEL_' time_stamp];
save(filename,'matlabbatch');                         

% Execute matlabbatch
if execDartel == 1
    spm_jobman('run',matlabbatch);
    disp(['DARTEL COMPLETED: ' datestr(now,'yyyymmdd_HHMM')]);
end
diary off