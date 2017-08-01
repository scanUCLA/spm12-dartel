function [status, errorMsg, allfuncs, allt1, allmt1, allrc1, allrc2, allu_rc1, allc1, allc2, allc3] =...
    run_dartel_mprageMBW(subNam, owd, codeDir, batchDir, runID, funcID, mbwdirID,...
    mpragedirID, fourDnii, execPreDartel)

% Last revision: 31 July 2017 - Kevin Tan
%% Parameters
funcFormat=2;       % format of your raw functional images (1=img/hdr, 2=4D nii)

%% Make Matlabbatch
spm('defaults','fmri'); spm_jobman('initcfg');

status = NaN;
errorMsg = '';
allfuncs = [];
allt1 = [];
allmt1 = [];
allrc1 = [];
allrc2 = [];
allu_rc1 = [];
allc1 = [];
allc2 = [];
allc3 = [];

try
    % Find run directories
    d = dir([owd '/' subNam '/raw/' runID]);
    nRuns = length(d);
    fprintf('Found %d runs\n', nRuns)
    tmpStr = sprintf('^%s.*\\.nii$', funcID);
    
    % Fund functional images per run
    volsFPr = {};
    for r = 1:nRuns
        vols = {};
        volsr = {};
        [vols, ~] = spm_select('ExtList', [d(r).folder '/' d(r).name], tmpStr, Inf);
        vols = cellstr(strcat(vols));
        for v = 1:length(vols)
            chars = length(vols{v})-2;
            volsFP{r}{v,1} = [d(r).folder '/' d(r).name '/' vols{v}];
            volsr{v,1} = [d(r).folder '/' d(r).name '/r' vols{v}(1:chars)];
        end 
        % Realigned images
        if fourDnii
            volsFPr = [volsFPr;volsr{1}];
        else
            volsFPr = [volsFPr;volsr];
        end  
    end
    
    % find the mbw folder
    d = dir([owd '/' subNam '/raw/' mbwdirID]);
    mbwdir = [d(1).folder filesep d(1).name];
    
    % find the mprage folder
    d = dir([owd '/' subNam '/raw/' mpragedirID]);
    mprdir = [d(1).folder filesep d(1).name];
    
    % get the images
    d = dir([mbwdir '/' mbwdirID '.nii']);
    mbw_name = d.name; clear d
    d = dir([mprdir '/' mpragedirID '.nii']);
    mprage_name = d.name; clear d
    mbw = [mbwdir filesep mbw_name];
    allt1 = [mprdir filesep mprage_name];
    fprintf('MBW is: %s\n',mbw)
    fprintf('MPRAGE is: %s\n\n',allt1)
    
    % for DARTEL
    allfuncs = volsFPr;
    allmt1 = [mprdir filesep 'm' mprage_name];
    allrc1 = [mprdir filesep 'rc1' mprage_name(1:end-4) '.nii'];
    allrc2 = [mprdir filesep 'rc2' mprage_name(1:end-4) '.nii'];
    allu_rc1 = [mprdir filesep 'u_rc1' mprage_name(1:end-4) '_Template.nii'];
    allc1 = [mprdir filesep 'c1' mprage_name];
    allc2 = [mprdir filesep 'c2' mprage_name];
    allc3 = [mprdir filesep 'c3' mprage_name];
    
    % =====================================
    % Begin building MATLABBATCH
    % =====================================
    
    matlabbatch{1}.cfg_basicio.cfg_cd.dir = cellstr(strcat([owd '/' subNam],filesep,'notes'));
    
    % Realignment of functionals
    % -------------------------------------------------
    for i = 1:nRuns
        matlabbatch{2}.spm.spatial.realign.estwrite.data{i} = volsFP{i};
    end
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;         % higher quality
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;               % default is 4
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;              % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0;               % 0=realign to first, to 1=realign to mean for
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 4;            % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];        % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = {};           % don't weight
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which  = [2 1];        % reslice all and mean
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;            % default
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap   = [0 0 0];      % no wrap (default)
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask   = 1;            % enable masking (default)
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Co-register MBW to MEAN FUNCTIONAL
    % -------------------------------------------------
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{3}.spm.spatial.coreg.estimate.source = cellstr(mbw);
    matlabbatch{3}.spm.spatial.coreg.estimate.other{1} = '';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % Co-register MPRAGE to MBW
    % -------------------------------------------------
    matlabbatch{4}.spm.spatial.coreg.estimate.ref = cellstr(mbw);
    matlabbatch{4}.spm.spatial.coreg.estimate.source = cellstr(allt1);
    matlabbatch{4}.spm.spatial.coreg.estimate.other{1} = '';
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
catch
    status = 0;
    errorMsg = 'Error making matlabbatch';
    disp([errorMsg ' for ' subNam]);
    cd(codeDir);
    return
end

%% Save matlabbatch
try
    time_stamp = datestr(now, 'yyyymmdd_HHMM');
    filename = [batchDir '/preDARTEL_MBWmprage_' subNam '_' time_stamp];
    save(filename, 'matlabbatch');
catch
    status = 0;
    errorMsg = 'Error saving matlabbatch';
    disp([errorMsg ' for ' subNam]);
    cd(codeDir);
    return
end
%% Run Matlabbatch
try
    if execPreDartel == 1
        spm_jobman('run',matlabbatch);
    end
    status = 1;
    cd(codeDir);
catch
    status = 0;
    errorMsg = 'Error running matlabbatch';
    disp([errorMsg ' for ' subNam]);
    cd(codeDir);
    return
end
end
