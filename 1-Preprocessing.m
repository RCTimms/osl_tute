clear all;
clc;
%% Set-up OSL and FSL
cd('/home/disk4/laurenza/osl-new/osl-core') % change this to reflect where your version of osl is downloaded
osl_startup;

%% Initialise file settings 
% Directory where the raw data is stored
rawdir = '/home/disk4/laurenza/Preprocessing_Workshop/Raw_Data/'; % change this to point to the raw data folder 
% Directory where the maxfiltered data is stored
datadir = '/home/disk4/laurenza/Preprocessing_Workshop/Maxfiltered_Data/'; % change this to point to the maxfiltered data folder

% Create an object including file path and name of raw data files
filePattern = fullfile(rawdir, '*.fif');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
  baseFileName = fifFiles(f).name;
  rawfif_files(f,:) = cellstr(fullfile(rawdir, baseFileName));
  clear baseFileName fullFileName
end
clear filePattern fifFiles

% Create an object including file path and name of maxfiltered data files
filePattern = fullfile(datadir, '*tsss.fif');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
  baseFileName = fifFiles(f).name;
  maxfif_files(f,:) = cellstr(fullfile(datadir, baseFileName));
  clear baseFileName fullFileName
end
clear filePattern fifFiles

%% Import raw fif to view data before maxfilter 

% Import the raw data (creates a .dat and .mat file in rawdir)
for subnum = 1:length(rawfif_files)
    D = osl_import(rawfif_files{subnum}); 
    clear D
end

filePattern = fullfile(rawdir, '*.mat');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
  baseFileName = fifFiles(f).name;
  rawmat_files(f,:) = cellstr(fullfile(rawdir, baseFileName));
  clear baseFileName fullFileName
end
clear filePattern fifFiles

for subnum = 1:length(rawmat_files)
    D = spm_eeg_load(rawmat_files{subnum});
    D = oslview(D)
    waitfor(gca)
    clear D
end

%% View maxfiltered data  

% Import the maxfiltered data (creates a .dat and .mat file in datadir)
for subnum = 1:length(maxfif_files)'
    D = osl_import(maxfif_files{jsubnum});
    clear D
end

filePattern = fullfile(datadir, '*.mat');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
  baseFileName = fifFiles(f).name;
  maxmat_files(f,:) = cellstr(fullfile(datadir, baseFileName));
  clear baseFileName fullFileName
end
clear filePattern fifFiles

% View the data to ensure maxfilter does something sensible and the data is relatively clean. 
% Also sanity check the spectra.

for subnum = 1:length(maxmat_files)
    D = spm_eeg_load(maxmat_files{subnum});
    D = oslview(D)
    waitfor(gca)
    osl_quick_spectra(D) % default is MEGPLANARS
    waitfor(gca)
    clear D
end

%% Inspect HPI fits for each file

filePattern = fullfile(datadir, '*hpi.log'); % may also be .pos file
logFiles = dir(filePattern);
for f = 1:length(logFiles)
  baseFileName = logFiles(f).name;
  fullFileName = fullfile(datadir, baseFileName);
  hpi_temp = osl_headpos(fullFileName);
  hpi(f,:) = hpi_temp;
  hpi_fits(f,:) = mean(hpi(f,:).hpi_fit); % values should be very high ie. > .99
  clear hpi_temp baseFileName fullFileName
end

% Add subplot function here for hpi_fit, translation and rotation shown
% over time? (subplot is currently in osl-core/+reports/headpos.m)

%% Co-registration

% In spm_sss we run the coregistration. This step does not interact with the MEEG data themselves
% so we keep it as a separate stage.
% If true, |use_existing| will prevent rerunning the coregstration on any files which already contain a coregistration

use_existing = true;

outdir = '/home/disk4/laurenza/Preprocessing_Workshop/Coregistration';    
mkdir(outdir)

% Ryan, I don't have any structural files for this data, will see if I can use COVID-MEG data!
load('/home/disk4/laurenza/Preprocessing_Workshop/StructuralMRI/StructuralFiles.mat')

% Main loop through subjects and sessions
for subnum = 1:length(maxmat_files)
    
    % Load data in from spm_sss, note that any changes in this loop are saved into the same file.
    D = spm_eeg_load(maxmat_files{subnum});
       
    % Coregistration is carried out using a call to osl_headmodel.
    coreg_settings = struct;
    coreg_settings.D = D.fullfile;
    coreg_settings.mri = structural_files(subnum);
    coreg_settings.useheadshape = true;
    coreg_settings.forward_meg = 'Single Shell';
    coreg_settings.use_rhino = true;
    coreg_settings.fid.label.nasion='Nasion';
    coreg_settings.fid.label.lpa='LPA';
    coreg_settings.fid.label.rpa='RPA';
    D = osl_headmodel(coreg_settings);
    
    % Next we generate and save a summary image so we can check each
    % coregistration has been carried out sucessfully.
    h = report.coreg(D);
    report.save_figs(h,outdir,D.fname);
    close(h);
    
    clear h D
end

%% Preprocessing
% Carried out on maxfiltered data files, imported to matlab as SPM objects 

% Downsampling - Downsample from 1000Hz to 250Hz
% Filtering - We apply a single passband filter to isolate data between 1 and 45Hz
% Bad Segment Detection - An automatic algorithm is used to identify noisy data segments which are removed from subsequent analysis
% Independant Components Analysis - ICA is used to identify artefactual componets by correlation with the EOG and ECG, these are removed from subsequent analysis
% Sensor Normalisation - The Magnetometers and Gradiometers within each dataset are normalised to make their variances comparable.
% Beamforming - An LCMV Beamformer is used to project the sensor data into an 8mm grid in source space
% Parcellation - The source space data is reduced into a network of parcels defined by a NIFTI file

% The downsampling and filtering overwrite the old data file, but the ICA,
% sensornormalisation (and beamforming and parcellation for source space data)
% are applied by adding online montages to the SPM object.

indir = datadir;
outdir = '/home/disk4/laurenza/Preprocessing_Workshop/Preprocessed_Data';
mkdir(outdir)
wd = '/home/disk4/laurenza/Preprocessing_Workshop/Temp_Dir';
mkdir(wd)

% Specify parcellation file for source space
p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );

sessions_to_check = []; % This will catch sessions to double check, although double checking of all data is recommended

% Note for Ryan - maybe best not to do the below in a loop, and go through
% one session in detail? 

% Main loop through files within the study object
for subnum = 1:length(maxmat_files)

	% Load in MEEG data as an SPM object
    D = spm_eeg_load(maxmat_files{subnum});

	% Downsample and copy, or just copy
	if D.fsample > 250
		D = spm_eeg_downsample(struct('D',D,'fsample_new',250,'prefix',[wd '/'])); % Note - downsampling cannot be done in-place using prefix='', it just fails
	else
		D = D.copy(getfullpath(fullfile(wd,D.fname))); % Copy into working directory
    end
    
	% Apply a 1-45Hz passband filter - could also be a low-pass or
	% high-pass filter
    D = osl_filter(D,[1 45],'prefix','');
    
    D = oslview(D); % view the data after downsampling/filtering to see effect
    waitfor(gca)
    osl_quick_spectra(D); % note effect of filtering on spectra 
    waitfor(gca)

	% Apply automatric bad segment detection
	D = osl_detect_artefacts(D,'badchannels',false,'badtimes',true); 

    D = oslview(D); % See where bad segments are marked in OSLview
    waitfor(gca)
    
	% Run ICA artefact detection. This will automatically reject components
	% which havoe correlations larger than .5 with either of the artefact
	% channels.
	D = osl_africa(D,'used_maxfilter',true);

    % Check warnings in the command line to see if no ICs are detected 
    % Check the D object now - you will see it has a montage applied.
    D
    
	% Though the automatic correlations generally work well, we should be
	% careful to check for unsusual datasets and possibly manually correct the
	% automatic assessment. This is particularly important for relatively noisy
	% data or when analysing a new dataset for the first time.
	
	% Here we will manally inspect the ICA component rejections for any dataset meeting the following criteria
	% # More than 4 ICs rejected
	% # Zero ICs rejected
	% # No component rejected due to correlation with EOG
	% # No component rejected due to correlation with ECG
    if isempty(D.ica.bad_components) || length(D.ica.bad_components) > 4
        disp('%s components rejected, recommend checking session', length(D.ica.bad_components));
        sessions_to_check = cat(1,sessions_to_check,subnum);
    elseif max(D.ica.metrics.corr_chan_308_EOG.value) < .5 || ...
       max(D.ica.metrics.corr_chan_307_ECG.value) < .5
        disp('no candidate components for either HEOG or ECG, recommend checking session');
        sessions_to_check = cat(1,sessions_to_check,subnum);
    end
    
    % This is where manual Africa would go
    D = D.montage('remove',1:D.montage('getnumber')); % This removes the montage from the automatic AFRICA
    D = osl_africa(D,'do_ident','manual'); % Montage with manual AFRICA will be saved

    % Normalise sensor types
    S = [];
    S.D = D;
    S.modalities = {'MEGMAG','MEGPLANAR'};
    S.do_plots = 0;
    S.samples2use = good_samples(D,D.indchantype(S.modalities,'GOOD'));
    S.trials = 1;
    S.pca_dim = 99;
    S.force_pca_dim = 0;
    S.normalise_method = 'min_eig';
    D = normalise_sensor_data( S );
   
    % Check the D object now - you will see it has 2 montages applied, the
    % latest being the sensor normalised data
    D
    
    % This is where the pipeline would end for data analysed in sensor space 
    
% 	% Run LCMV Beamformer
% 	D = osl_inverse_model(D,p.template_coordinates,'pca_order',50);
% 
% 	% Do parcellation
% 	D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
% 	
	% Save out
    D.save
	D = D.montage('switch',0);
	D.copy(fullfile(outdir,D.fname));
    
    clear D
end