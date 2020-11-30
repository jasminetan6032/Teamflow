function connectivity_results = get_connectivity(sub) 
%This script takes EEGlab data from two partners and extracts inter-brain connectivity measures 
%input variable "sub" should be initials for two participants (take from
%participant list)eg. ['AL' ; 'TB']

%indexing 2 participants in data directory
whereis = 'E:\Jasmine\MSc_in_Music_Mind_and_Brain\Research Project\EEG\Analysed EEG\Analysed Data (from lab computer)\InterpEpoch80\PF_FLOW_80_select\';
findset = dir([whereis '\*.set']); % find .set files
findset = findset(~[findset.isdir]');
partnames = struct2cell(findset);
partnames = transpose(partnames(1,:));
index = cell2mat(cellfun(@(x) any(strncmp(x,sub(1,:),2) || strncmp(x,sub(2,:),2)), partnames,'UniformOutput',false));

% hello

%gives error if participants are not found 
if nargin <1 || isnumeric(sub) 
    disp ('Error: Please input two initials')
elseif sum(index)~= 2
    disp ('Participant does not exist')
end

%loads EEGlab 
filenames = findset(index);
addpath('D:\EEG_analyses\eeglab versions\eeglab_current\eeglab2019_1');
load('D:\EEG_analyses\EEG\2020 Experiments\Group flow\customchanlocs.mat');

%load both datasets
EEG1 = pop_loadset('filename',filenames(1).name,'filepath',whereis);
EEG1 = pop_select( EEG1,'channel',{'Fp1' 'AF7' 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'AFz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'});

EEG2 = pop_loadset('filename',filenames(2).name,'filepath',whereis);
EEG2 = pop_select( EEG2,'channel',{'Fp1' 'AF7' 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'AFz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'});

%combine both EEG
combine_EEG = cat(1, EEG1.data,EEG2.data);

%replace EEG1 with combined data
EEG1.data = combine_EEG;

%replace EEG1 chanlocs with customchanlocs for 2 electrode sets
EEG1.chanlocs = customchanlocs;

%convert to fieldtrip to do connectivity analysis
EEG_combine_ft = eeglab2fieldtrip(EEG,'preprocessing','none');
addpath(genpath('E:\Jasmine\MSc_in_Music_Mind_and_Brain\fieldtrip-20140401'));

%Time-frequency analysis 
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi       = 1:1:70;
cfg.tapsmofrq = cfg.foi*0.05;
freqdata = ft_freqanalysis(cfg,EEG_combine_ft);

%Connectivity analysis
cfg = [];
cfg.method = 'psi'; %this can be changed for other connectivity measures eg. wpli_debiased, plv
cfg.bandwidth = 1; % change here according to the frequency resolution
connectres_combine(part_i) = ft_connectivityanalysis(cfg,freqdata);

%get part of matrix that corresponds to interbrain connectivity (note: from
%participant 1 to participant 2)
connectivity_results = connectres_combine.psispctrm(1:64,65:128,:);

rmpath(genpath('E:\Jasmine\MSc_in_Music_Mind_and_Brain\fieldtrip-20140401'));


%% 3. Plot the head-in-head using the heads-bio function
% 
% % this first bit is to load the needed locations
% load elec_biosemi64; %load some eeg 3D-channels locations; 31 channels 
% clear para; para.rot=180; %setting this rotates locations by 180 degrees
% elec = elec_biosemi64.pnt(1:64,:);
% locs_2D=mk_sensors_plane(elec); 
% 
% % define the y limits as
% cfg.ylim = [-0.05 0.05];
% cfg.freqband = [5 7.5];
% % feed that in the heads_bio function
% heads_bio(gavgCON,cfg.freqband,cfg.ylim,locs_2D,1)

end
