%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Novelty Oddball Parametric Modulation Constrats  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for setting up a contrast for the parametric modulation analysis.

%adjust path were all the needed data is stored
studyPath = '/Users/danieljanko/Desktop/Projects/Oddball/oddball/'; % where your data are

% Where to save the model specification
parametric_modulation = '/Param_modul';

% loading the table with IDs
subjects = readtable('/Users/danieljanko/Desktop/Projects/Oddball/oddball/participants.tsv', 'FileType', 'text');
n_sub = height(subjects);

for subj = 1:n_sub
    ID = ['sub-' char(subjects{subj,1})];
    directory = fullfile(studyPath, 'derivatives', ID);
    spm_mat = {strcat(directory, parametric_modulation, '/SPM.mat')};
    matlabbatch{1}.spm.stats.con.spmmat = spm_mat;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'All';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1 1 1 1 1 1 1	1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0	0	0	0	0 0]; 
    spm_jobman('run', matlabbatch(1))
end 


