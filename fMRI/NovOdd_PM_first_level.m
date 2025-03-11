%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         First Level Parametric Modulation         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up a first level model for the parametric modulation analysis. 

% Clear and close all 
clear;
clc;

% Initialise SPM
spm('defaults','fmri');  
spm_jobman('initcfg');

%adjust path were all the needed data is stored
studyPath = '/Users/danieljanko/Desktop/Projects/Oddball/oddball/'; % where your data are

% Where to save the model specification
parametric_modulation = '/Param_modul';

% loading the table with IDs
subjects = readtable('/Users/danieljanko/Desktop/Projects/Oddball/oddball/participants.tsv', 'FileType', 'text');
n_sub = height(subjects);
TR = 2;

for subj = 1:n_sub
    ID = ['sub-' char(subjects{subj,1})];
    directory = fullfile(studyPath, 'derivatives', ID);
    mkdir(sprintf('%s/Param_modul', directory))
    save = {strcat(directory, parametric_modulation)}; % where is the first level folder saved
    matlabbatch{1}.spm.stats.fmri_spec.dir = save;
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    logdir = directory;
    onset = strcat(logdir, '/', ID, '_task-oddball_acq-2mmMB2S2_events.tsv');
    onsetfile = readtable(onset, 'FileType', 'text', 'Delimiter', 'tab');
    col = find(onsetfile{:, 1} == 3 | onsetfile{:, 1} == 4);
    ons = table2array(onsetfile(col, 2));
    for c = 1:length(col)
        param_modulator = round(str2double(onsetfile{col(c), "decay"}), 3);
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).name = 'Onset Oddball' ;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).onset = ons(c); 
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).orth = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).pmod(1).name = ['Oddball' num2str(c)];  % Parametric modulator name (e.g., reaction time)
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).pmod(1).param = param_modulator;  % Latent variable for pre-training
    end 
    fun_file = cellstr(spm_select('ExtFPList', directory, '^swar.*\.nii$', 1:365));
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(fun_file);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    rp = spm_select('FPList', directory, '^rp_.*\.txt$');
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(rp);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; 
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    spm_jobman ('run', matlabbatch(1));

    %%%%%%% ESTIMATE %%%%%%%%

    % Beta estimation
    spm_mat = {strcat(directory, parametric_modulation, '/SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.spmmat = spm_mat;
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
    % Run batch n2 
    spm_jobman ('run', matlabbatch(2));
end 