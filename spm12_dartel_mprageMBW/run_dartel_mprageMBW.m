function [status, errorMsg, allfuncs, allt1, allmt1, allrc1, allrc2, allu_rc1] =...
    run_dartel_mprageMBW(subNam, owd, codeDir, batchDir, runID, funcID, mbwdirID,...
    mpragedirID, fourDnii, execPreDartel)

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

try
    cd(owd)
    swd = sprintf('%s/%s',owd,subNam);
    fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd raw
    base_dir = pwd;
    
    % Find run directories
    d=dir(runID);
    run_names = {d.name};
    numruns=length(run_names);
    fprintf('Found %d runs\n',numruns)
    
    % Find functional images for run(s)
    %-----------------------------------------------------------------%
    load_dir = {};
    raw_func_filenames = {};
    allfiles_orig = {};
    allfiles_norm = {};
    for i = 1:numruns
        load_dir{i} = fullfile(base_dir,run_names{i});
        if funcFormat==1
            tmpString=sprintf('^%s.*\\.img$',funcID);
            [raw_func_filenames{i},dirs] = spm_select('List',load_dir{i},tmpString, inf);
            filenames_orig{i}=cellstr(strcat(load_dir{i},filesep,raw_func_filenames{i}));
            filenames_norm{i}=cellstr(strcat(load_dir{i},filesep,'w',raw_func_filenames{i}));
            allfiles_orig = [allfiles_orig; filenames_orig{i}];
        else
            tmpString=sprintf('^%s.*\\.nii$',funcID);
            [raw_func_filenames{i},dirs] = spm_select('ExtFPList',load_dir{i},tmpString, inf);
            filenames_orig{i}=cellstr(strcat(raw_func_filenames{i}));
            filenames_norm{i}=cellstr(strcat('w',raw_func_filenames{i}));
            allfiles_orig = [allfiles_orig; filenames_orig{i}];
        end
    end
    if funcFormat==1
        mean_func=cellstr(strcat(load_dir{1},filesep,'mean',raw_func_filenames{1}(1,:)));
    else
        [path name ext] = fileparts(allfiles_orig{1});
        mean_func=cellstr(strcat(path,filesep,'mean',name,'.nii'));
    end
    load_dir = fullfile(base_dir,run_names{i});
    
    % Find the anatomicals
    % -------------------------------------------------
    
    % find the mbw folder
    d=dir(mbwdirID);
    mbwdir = [base_dir filesep d(1).name];
    
    % find the mprage folder
    d=dir(mpragedirID);
    mprdir = [base_dir filesep d(1).name];
    
    % get the images
    cd(mbwdir); d = dir([mbwdirID '.nii']);
    mbw_name = d.name; clear d
    cd(mprdir); d = dir([mpragedirID '.nii']);
    mprage_name = d.name; clear d
    mbw = [mbwdir filesep mbw_name];
    allt1 = [mprdir filesep mprage_name];
    fprintf('MBW is: %s\n',mbw)
    fprintf('MPRAGE is: %s\n\n',allt1)
    [firstpart,lastpart] = strread(allt1,'%s %s','delimiter','.');
    
    % for DARTEL
    if fourDnii == 1
        allfuncs = filenames_orig;
    else
        allfuncs = allfiles_orig;
    end
    allmt1 = [mprdir filesep 'm' mprage_name];
    allrc1 = [mprdir filesep 'rc1' mprage_name(1:end-4) '.nii'];
    allrc2 = [mprdir filesep 'rc2' mprage_name(1:end-4) '.nii'];
    allu_rc1 = [mprdir filesep 'u_rc1' mprage_name(1:end-4) '_Template.nii'];
    
    % =====================================
    % Begin building MATLABBATCH
    % =====================================
    
    matlabbatch{1}.cfg_basicio.cfg_cd.dir = cellstr(strcat(swd,filesep,'notes'));
    
    % Realignment of functionals
    % -------------------------------------------------
    for i = 1:numruns
        matlabbatch{2}.spm.spatial.realign.estwrite.data{i} = filenames_orig{i};
    end
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;         % higher quality
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;               % default is 4
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;              % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;               % changed from 0 (=realign to first) to 1 (realign to mean) for
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 4;            % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];        % default
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = {};           % don't weight
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which  = [0 1];        % create mean image only when reslicing
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;            % default
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap   = [0 0 0];      % no wrap (default)
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask   = 1;            % enable masking (default)
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Co-register MBW to MEAN FUNCTIONAL
    % -------------------------------------------------
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = mean_func;
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
