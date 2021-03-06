%% House keeping!
clear all; clc; close all;
%% Set-up OSL and FSL paths
cd('/home/disk4/laurenza/osl-new/osl-core') % change this to reflect where your version of osl is downloaded
osl_startup;

%% Initialise file settings
% This cell just sets up some objects pointing to where all our data is
% stored. This can be really useful if you have many subjects' data

% Directory where the raw data is stored. You can change this to point to
% your raw data folder!

rawdir = '/Volumes/TASER/Lauren_osl_course/henson_sub01/raw';


% Directory where your maxfiltered data is stored (if you have already run
% it)
datadir = '/Volumes/TASER/Lauren_osl_course/henson_sub01/maxfilter_output';

% Create an object including file path and name of raw data files.
filePattern = fullfile(rawdir, '*.fif');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
    baseFileName = fifFiles(f).name;
    rawfif_files(f,:) = cellstr(fullfile(rawdir, baseFileName));
    clear baseFileName fullFileName
end
clear filePattern fifFiles

% Create an object including file path and name of maxfiltered data files.
% This sometimes has the extension *tsss, sometimes it's just *sss
filePattern = fullfile(datadir, '*sss.fif');
fifFiles = dir(filePattern);
for f = 1:length(fifFiles)
    baseFileName = fifFiles(f).name;
    maxfif_files(f,:) = cellstr(fullfile(datadir, baseFileName));
    clear baseFileName fullFileName
end
clear filePattern fifFiles

%% Import raw fif to view data before maxfilter
% Let's see what our data looks like before we do anything to it (brace
% yourself!). We'll use osl_import to read in the data, followed by oslview
% to visualise it.

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
% Now let's see what effect the MaxFilter operation has on the data.

% Import the maxfiltered data (creates a .dat and .mat file in datadir)
for subnum = 1:length(maxfif_files)'
    D = osl_import(maxfif_files{subnum});
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

% View the data to ensure maxfilter does something sensible and the data is
% relatively clean. Also sanity check the spectra (the frequency
% representation of the data). What do you notice about this?

for subnum = 1:length(maxmat_files)
    D = spm_eeg_load(maxmat_files{subnum});
    D = oslview(D)
    waitfor(gca)
    
    % We have two options for visualising the spectra. To make things
    % efficient here we'll call our own pwelch function on just one
    % channel, but once we've downsampled we'll use osl_quick_spectra(D)
    data=D(:,:,:);
    figure;pwelch(data(3,1:1e4),[],[],[],D.fsample);
    xlabel('Frequency (Hz)');ylabel('PSD');
    title(['PSD for channel ',D.chanlabels(3)]);xlim([1 70]);
    set(gca,'FontSize',16)
    
    waitfor(gca)
    clear D
end

% Note that it might be a good idea to save-out the "max_mat" files
% somewhere to disk, in case you want to pause your analysis or only re-run
% part of the pipeline.
%% Explore the MaxFiltered D object
% Let's load in our MaxFiltered data. What do you notice about the D
% object? What can you do with it? Have a play!
D = spm_eeg_load(maxmat_files{1})


%% Renaming of peripherary channels.
% This shouldn't be an issue when working on Oxford data. Here we have a
% dataset collected at Cambridge and so we have to relabel some of the EEG
% channels to their correct artefact names. We have two ocular and one
% heartbeat recording.
D = D.chantype(find(strcmp(D.chanlabels,'EEG062')),'EOG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG062')),'VEOG');

D = D.chantype(find(strcmp(D.chanlabels,'EEG061')),'EOG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG061')),'HEOG');

D = D.chantype(find(strcmp(D.chanlabels,'EEG063')),'ECG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG063')),'ECG');

D.save();

%% Inspect HPI fits for each file
% MaxFilter does a lot of useful operations which are beyond the scope of
% this current tutorial. However, it also provides a useful output file
% which we can run some sanity checks on to ensure that it has run
% successfully. We can also track the subjects' head position throughout an
% experimental recording - we might want to discard any subjects who fidget
% too much!

% The following line should work on data collected in Oxford: filePattern =
% fullfile(datadir, '*hpi.log'); % may also be .pos file logFiles =
% dir(filePattern); for f = 1:length(logFiles)
%   baseFileName = logFiles(f).name; fullFileName = fullfile(datadir,
%   baseFileName); hpi_temp = osl_headpos(fullFileName,1); hpi(f,:) =
%   hpi_temp; hpi_fits(f,:) = mean(hpi(f,:).hpi_fit); clear hpi_temp
%   baseFileName fullFileName
% end


% The naming convention in Cambridge is slightly different, so we'll just
% take that as given here:
hpos_name = '/Volumes/TASER/Lauren_osl_course/henson_sub01/run_01_headpos.txt';
hpi_temp = osl_headpos(hpos_name,1);
hpi(f,:) = hpi_temp;
hpi_fits(f,:) = mean(hpi(f,:).hpi_fit); % values should be very high ie. > .99


%% Co-registration

% In spm_sss we run the coregistration. This step does not interact with
% the M/EEG data themselves so we keep it as a separate stage. If true,
% |use_existing| will prevent rerunning the coregstration on any files
% which already contain a coregistration. This step can take a little
% while. Note that RHINO/SPM will output quite a few additional MRI-like
% files in the same directory as the original structural file.

use_existing = true;

% outdir = '/home/disk4/laurenza/Preprocessing_Workshop/Coregistration';
% mkdir(outdir)

structural_files={'/Volumes/TASER/Lauren_osl_course/henson_sub01/mprage.nii'}

% Main loop through subjects and sessions
for subnum = 1:length(maxmat_files)
    
    % Load data in from spm_sss, note that any changes in this loop are
    % saved into the same file.
    D = spm_eeg_load(maxmat_files{subnum});
    
    % Coregistration is carried out using a call to osl_headmodel.
    coreg_settings = struct;
    coreg_settings.D = D.fullfile;
    coreg_settings.mri = structural_files(subnum);
    coreg_settings.useheadshape = false;
    coreg_settings.forward_meg = 'Single Shell';
    coreg_settings.use_rhino = true;
    coreg_settings.fid.label.nasion='Nasion';
    coreg_settings.fid.label.lpa='LPA';
    coreg_settings.fid.label.rpa='RPA';
    D = osl_headmodel(coreg_settings);
    
    % LOOK AT YOUR DATA! Call out the RAC.
    osl_RAC(D);
    clear h D
end

%% Preprocessing
% Carried out on maxfiltered data files, imported to matlab as SPM D
% objects. We have several key pre-processing steps to carry out.

% 1) Downsampling - Downsample the data from 1000Hz to 250Hz. If you're
% wanting to inspect gamma-band activity then this will likely have to be
% higher than 250Hz (see discussion).

% 2) Filtering - Here we apply a single passband filter to isolate data
% between 1 and 45Hz. This is likely to change depending on your hypothesis
% (e.g. slow rhythms, motor beta band activity, etc.).

% 3) Notch filter - Filters aren't always perfect, so it's best to err on
% the side of caution and throw in a secondary 50Hz filter to remove any
% residual line noise.

% 4) Bad Segment Detection - An automatic algorithm is used to identify
% noisy data segments which are removed from subsequent analysis.

% 5) Independant Components Analysis - ICA is used to identify artefactual
% componets by correlation with the EOG and ECG, these are removed from
% subsequent analysis. This often requires user input.

% 6) Sensor Normalisation - The Magnetometers and Gradiometers within each
% dataset are normalised to make their variances comparable (if you have
% Elekta data).

% 7) Beamforming - An LCMV Beamformer is used to project the sensor data
% onto an 8mm grid in source space. Other grid sizes are available!

% 8) Parcellation - The source space data is reduced into a network of
% parcels defined by a NIFTI file. This can be thought of like an atlas.

% The downsampling and filtering overwrite the old data file, but the ICA,
% sensornormalisation (and beamforming and parcellation for source space
% data) are applied by adding online montages to the SPM object.

indir = datadir;
outdir = '/Volumes/TASER/Lauren_osl_course/henson_sub01/Preprocessed_data';
mkdir(outdir)
wd = '/Volumes/TASER/Lauren_osl_course/henson_sub01/tmp_dir';
mkdir(wd)

% Specify parcellation file for source space
p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );

sessions_to_check = []; % This will catch sessions to double check, although double checking of all data is recommended


% Main loop through files within the study object
for subnum = 1:length(maxmat_files)
    
    % Load in MEEG data as an SPM object
    D = spm_eeg_load(maxmat_files{subnum});
    
    % 1) Downsample and copy, or just copy D object
    if D.fsample > 250
        D = spm_eeg_downsample(struct('D',D,'fsample_new',250,'prefix',[wd '/'])); % Note - downsampling cannot be done in-place using prefix='', it just fails
    else
        D = D.copy(getfullpath(fullfile(wd,D.fname))); % Copy into working directory
    end
    
    % 2) Apply a 1-45Hz passband filter - could also be a low-pass or
    % high-pass filter
    D = osl_filter(D,[1 45],'prefix','');
    D = oslview(D); % view the data after downsampling/filtering to see effect
    osl_quick_spectra(D); % note effect of filtering on spectra
    drawnow
    
    % 3) Also apply a 48-52Hz notch filter
    D = osl_filter(D,-[48 52],'prefix','');
    D = oslview(D); % view the data after downsampling/filtering to see effect
    osl_quick_spectra(D); % note effect of filtering on spectra
    drawnow
    
    % 4) Apply automatric bad segment detection
    D = osl_detect_artefacts(D,'badchannels',false,'badtimes',true);
    
    D = oslview(D); % See where bad segments are marked in OSLview
    waitfor(gca)
    
    % 5) Run ICA artefact detection. This will automatically reject
    % components which have correlations larger than .5 with either of the
    % artefact channels.
    D = osl_africa(D,'used_maxfilter',true);
    
    % Check warnings in the command line to see if no ICs are detected
    % Check the D object now - you will see it has a montage applied.
    D
    
    % Though the automatic correlations generally work well, we should be
    % careful to check for unsusual datasets and possibly manually correct
    % the automatic assessment. This is particularly important for
    % relatively noisy data or when analysing a new dataset for the first
    % time. In short, we would **always** recommend **LOOKING AT YOUR
    % DATA** and take nothing for granted!
    
    % Note that the manual inspection of artefactual components can be
    % quite labour intensive. You might want to construct your for loop
    % over subjects in such a way that you only process 5 or 10 subjects at
    % a time so that you can have a break. We'd also recommend buying some
    % caffeine (coffee), fruit (satsumas) and snacks (pringles) to keep you
    % focussed!
    
    % Here we will manally inspect the ICA component rejections for any
    % dataset meeting the following criteria:
    %   1) More than 4 ICs rejected 2) Zero ICs rejected 3) No component
    %   rejected due to correlation with EOG(s) 4) No component rejected
    %   due to correlation with ECG
    
    % Note that we have had to change the channel numbers for the following
    % line to work as we have Cambridge data. If you were analysing data
    % collected in Oxford, uncomment out the following lines:
    
    %     if isempty(D.ica.bad_components) || length(D.ica.bad_components)
    %     > 4
    % fprintf('%d components rejected, recommend checking session\n',
    % length(D.ica.bad_components))
    %         sessions_to_check = cat(1,sessions_to_check,subnum);
    %     elseif max(D.ica.metrics.corr_chan_308_EOG.value) < .5 || ...
    %        max(D.ica.metrics.corr_chan_307_ECG.value) < .5
    %         disp('no candidate components for either HEOG or ECG,
    %         recommend checking session'); sessions_to_check =
    %         cat(1,sessions_to_check,subnum);
    %     end
    
    
    if isempty(D.ica.bad_components) || length(D.ica.bad_components) > 4
        fprintf('%d components rejected, recommend checking session\n', length(D.ica.bad_components))
        sessions_to_check = cat(1,sessions_to_check,subnum);
    elseif max(D.ica.metrics.corr_chan_367_EOG.value) < .5 || ...
            max(D.ica.metrics.corr_chan_368_EOG.value) < .5 || ...
            max(D.ica.metrics.corr_chan_369_ECG.value) < .5 || ...
            disp('no candidate components for either HEOG or ECG, recommend checking session');
        sessions_to_check = cat(1,sessions_to_check,subnum);
    end
    
    % This is where manual Africa would go
    D = D.montage('remove',1:D.montage('getnumber')); % This removes the montage from the automatic AFRICA
    D = osl_africa(D,'do_ident','manual'); % Montage with manual AFRICA will be saved
    
    % 6) Normalise sensor types
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
    
    % This is where the pipeline would end for resting state data analysed
    % in sensor space
    
    % 7) Run LCMV Beamformer
    D = osl_inverse_model(D,p.template_coordinates,'pca_order',50);
    
    % 8) Do parcellation
    D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
    
    % Save out
    D.save
    D = D.montage('switch',0);
    D.copy(fullfile(outdir,D.fname));
    
end
%%
% To finish, one of the key advantages of beamforming to an MNI-aligned
% grid is that we can use a bunch of built in tools to visualise our source
% space data on a template brain. This makes making figures for publication
% nice and easy!

% Switch back to the beamformer source space montage
D=D.montage('switch',4);

% Saving a nifti is as simple as calling the following:
out_fname=p.savenii(D(:,1),'my_beautiful_figure.nii');

% Here the output, out_fname, is the name of our new nifti file. The
% input is a nvoxels x time matrix. Here we have just given the first time
% point of our source space data, but we can also provide multiple time
% points, e.g.:

out_fname_movie=p.savenii(D(:,1:10),'my_beautiful_dynamic_figure.nii');

% We can then visualise these outputs in fsleyes:
fsleyes({out_fname,out_fname_movie})

% There are many useful things you can do with a parcellation object and
% the fsleyes function used above. Have a look here
% https://ohba-analysis.github.io/osl-docs/matlab/osl_example_parcellation.html#!
% but also just have a play!


% Not covered today: epoching