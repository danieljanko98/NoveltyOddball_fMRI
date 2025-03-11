%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Second Level Parametric Modulation        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for the second level parametric modulation

clear;
clc;

% Initialise SPM
spm('defaults','fmri');  
spm_jobman('initcfg');

%adjust path were all the needed data is stored
studyPath = '/Users/danieljanko/Desktop/Projects/Oddball/oddball/derivatives'; % where your data are
subjects = readtable('/Users/danieljanko/Desktop/Projects/Oddball/oddball/participants.tsv', 'FileType', 'text');
subjects = string(subjects{:,1});
file_paths = cell(height(subjects), 1);
for i = 1:height(subjects)
    file_paths{i} = sprintf('/Users/danieljanko/Desktop/Projects/Oddball/oddball/derivatives/sub-%s/Param_modul/con_0001.nii,1', subjects{i});
end
% Where to save the model specification
parametric_modulation = '/Param_modul';

dir = fullfile(studyPath, 'PM_second_level');
mkdir(dir);

matlabbatch{1}.spm.stats.factorial_design.dir = {dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {file_paths};

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch(1))


mat = fullfile(studyPath, 'PM_second_level/SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.spmmat = {mat};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch(2))