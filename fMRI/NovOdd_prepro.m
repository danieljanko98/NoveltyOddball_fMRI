%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Novelty Oddball fMRI Preprocessing       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for preprocessing of fMRI data of the Novelty Oddball task. This
% script is built to work with BIDS structured data. The raw fodler
% includes .nii.gz data that are unziped and the new file is transported
% into the derivatives folder. Further preprocessing steps are carried out
% on the moved file inside the derivatives directory.


% Clear and close all 
clear;
clc;

% Initialise SPM
spm('defaults','fmri');  
spm_jobman('initcfg');

%adjust path were all the needed data is stored
studyPath = '/Users/danieljanko/Desktop/Projects/Oddball/oddball/'; % where your data are

% loading the table with IDs
subjects = readtable('/Users/danieljanko/Desktop/Projects/Oddball/oddball/participants.tsv', 'FileType', 'text');
n_sub = height(subjects);
TR = 2;

for subj = 1:n_sub
    ID = ['sub-' char(subjects{subj,1})];
    fun_dir = [studyPath, ID, '/func/'];
    gunzip(fullfile(fun_dir, '*.nii.gz'));
    file = spm_select('ExtFPList', fun_dir, '.*bold.nii', 1);
    nif = nifti(file);
    N = nif.dat.dim(4);

    derivatives = fullfile(studyPath, 'derivatives', ID);
    if ~exist(derivatives, 'dir')
        mkdir(derivatives);
    end
    movefile(fullfile(fun_dir, '*bold.nii'), derivatives);

    %% SPATIAL REALIGNMENT %% 
    fun_file = spm_select('ExtFPList', derivatives, '.*bold.nii', 1:N);
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(fun_file)};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r_';

    spm_jobman('run', matlabbatch(1))

    
    %% SLICE TIMING CORRECTION %%
    v = spm_select('ExtFPList', derivatives, '^r_.*.nii', 1:N);
    V = spm_vol(v);
    Dim = V.dim;
    TA = TR - (TR/Dim(1,3));
    fun_file2 = spm_select('ExtFPList', derivatives, '^r_.*bold.nii', 1:N);
    matlabbatch{2}.spm.temporal.st.scans = {cellstr(fun_file2)};
    matlabbatch{2}.spm.temporal.st.nslices = Dim(1,3);
    matlabbatch{2}.spm.temporal.st.tr = TR;
    matlabbatch{2}.spm.temporal.st.ta = TA;
    matlabbatch{2}.spm.temporal.st.so = [1:2:Dim(1,3) 2:2:Dim(1,3)];
    matlabbatch{2}.spm.temporal.st.refslice = 1;
    matlabbatch{2}.spm.temporal.st.prefix = 'a';

    spm_jobman('run', matlabbatch(2))

    %% COREGISTER %%
    anat_dir = [studyPath, ID, '/anat/'];
    gunzip(fullfile(anat_dir, '*.nii.gz'));
    movefile(fullfile(anat_dir, '*T1w.nii'), derivatives);
    anatomy = spm_select('ExtFPList', derivatives, 'T1w.*.nii', 1);
    mean = spm_select('ExtFPList', derivatives, 'mean.*.nii', 1);

    matlabbatch{3}.spm.spatial.coreg.estwrite.ref = {mean};
    matlabbatch{3}.spm.spatial.coreg.estwrite.source = {anatomy};
    matlabbatch{3}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    spm_jobman('run', matlabbatch(3))

    %% SEGMENTATION %%
    anatomy2 = spm_select('ExtFPList', derivatives, '^r.*T1w.nii', 1);
    matlabbatch{4}.spm.spatial.preproc.channel.vols = {anatomy2};
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,1'};
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,2'};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,3'};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,4'};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,5'};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {'/Users/danieljanko/Documents/MATLAB/spm12/tpm/TPM.nii,6'};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{4}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{4}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                  NaN NaN NaN];

    spm_jobman('run', matlabbatch(4))

    %% NORMALIZE %% 

    files = dir(fullfile(derivatives, '*'));
    regex = '^y_.*\.nii$';
    def = fullfile(derivatives, {files(~cellfun('isempty', regexp({files.name}, regex))).name});
    fun_file3 = spm_select('ExtFPList', derivatives, '^ar_.*bold.nii', 1:N);

    matlabbatch{5}.spm.spatial.normalise.write.subj.def = cellstr(def{1});
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(fun_file3);
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run', matlabbatch(5))


    %% SMOOTH %%
    w = ['war_', ID, '.*.nii'];
    wfimages = spm_select('ExtFPList', derivatives, w, 1:N);

    matlabbatch{6}.spm.spatial.smooth.data = cellstr(wfimages);
    matlabbatch{6}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{6}.spm.spatial.smooth.dtype = 0;
    matlabbatch{6}.spm.spatial.smooth.im = 0;
    matlabbatch{6}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run', matlabbatch(6))

end 