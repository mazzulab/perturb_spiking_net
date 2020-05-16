% Demo script for running spiking network simulations and analyses
% 
% by Luca Mazzucato 2020
%
% ----------------------------------------
% Please cite:
% L. Mazzucato, G. La Camera, A. Fontanini 
% Expectation-induced modulation of metastable activity underlies faster coding of sensory stimuli, 
% Nat. Neuro. 22, 787-796 (2019).
% ----------------------------------------
% This script runs simulations of LIF networks with excitatory (E) and
% inhibitory (I) spiking neurons, with E and I clusters  [from Wyrick et al. 2020]

% CHOOSE STIMULI (OR CREATE YOUR OWN)
% You have 3 built-in stimulus templates:
% 0) no stimuli: run the network during ongoing activity
% 1) 'US' stimulus onset at t=0s (stimulus is on until the end of the trial);
%   a linearly ramping external current delivered to 50% of E clusters (chosen randomly, 
%   half of the neurons in each stimulus-selective clusters receive stimulus) with a slope of 
%   gain=0.2*mu_ext, where mu_ext is the baseline external current (bias)
% 2) 'CSgauss' stimulus onset at t=-0.5s (stimulus is on until the end of
%   the trial); double exponential profile with rise and decay times
%   [0.5,1]s; the stimulus targets all E neurons. For each neuron, the
%   stimulus peak is drawn from a gaussian distribution with mean 0 and
%   standard deviation 0.2*mu_ext. This is a contextual modulation
%   reproducing the \delta var(E) perturbation in Wyrick et al., 2020, with
%   CV(E)=10%.
% CUSTOM OPTIONS:
% 1) edit your custom-made stimuli inside aux.create_params. Available
% options: 'US', 'CSgauss' from Wyrick et al., 2020
% 2) list which events to run in current trials in the cell array 'events'
% 3) if events={}, run trial with spontaneous activity only, without any stimulus
%------------------------
% NOTE:
% This demo runs the unperturbed (UT) and perturbed (ET) conditions, depending on the following settings:
% 1) with condition='UT' you simulate an unperturbed trial, with a sensory
% stimulus delivered at t=0s and no contextual perturbation (Fig. 2a left panel in the paper)
% 2) with condition='ET' you simulate an expected trial, with a sensory stimulus
% delivered at t=0, preceded by a contextual perturbation at t=-0.5s (Fig. 2a right panel in the paper).
%-----------------
% SELECT STIMULUS
%-----------------
% events={}; % ongoing activity
% events={'US'}; % stimulus evoked-activity targeting selective clusters
% events={'CSgauss'}; % contextual perturbation speeds up network dynamics
%------------------------
ClustersOption='EI';%
%------------------------
% LOAD PARAMETERS
%------------------------
paramsfile='params.mat'; % file where all network parameters are saved
aux.create_params_EI(paramsfile);events={'US','CSgauss'}; % anticipatory cue preceeds stimulu delivery
save(paramsfile,'events','-append');
savedir=fullfile('data'); if ~exist(savedir,'dir'); mkdir(savedir); end % setup directory for saving HMM data

%% RUN SIMULATION
ntrials=2; % number of trials
file_sim=fullfile(savedir,'results.mat');  % file where simulation results are saved
%---------------------------
% GENERATE SYNAPTIC WEIGHTS
%---------------------------
% J = N x N matrix of synaptic weights
% params = structure containing all network parameters
[J, params]=aux.fun_SynWeights_EI(paramsfile);
[stimulus_save, params]=aux.fun_stim(params); % STIMULUS
%------------------------
% SIMULATION
%------------------------
tic
firings=cell(1,ntrials); % cell array with all spike times in each trial
PlotData=cell(1,ntrials); % cell array with data for plotting
for iTrial=1:ntrials
    ParamsRun=params;
    ParamsRun.Ext=stimulus_save.Ext;
    ParamsRun.Stimulus=stimulus_save.Stimulus;
    ParamsRun.J=J;
    fprintf('--- Start SIM ...\n');
    [firings{iTrial}, PlotData{iTrial}]=aux.fun_LIF_SIM(ParamsRun);
end
% SAVE results
save(file_sim,'params','firings','PlotData','stimulus_save');
fprintf('\nDone. Simulation saved in %s\n',file_sim);
toc
%%
%------------------------
% PLOT EVENTS
%------------------------
iTrial=2; % pick which trial to plot
dataload=load(file_sim); % load simulation results
data=dataload.PlotData{iTrial}; % membrane potentials
firings=dataload.firings{iTrial}; % spikes
Params=dataload.params; % parameters
Params.Ext=dataload.stimulus_save.Ext; % external currents
Params.savedir=savedir;
aux.fun_PlotTrial(data,firings,Params);
% figure 1 - rasterplot of all neurons in trial
% figure 2 - time course of membrane potential and PSC traces for E and I representative neurons
% figure 3 - time course of firing rate in clusters
% figure 4 - time course of stimuli (with CSgauss stimulus, the cue profile should be multiplied by a factor drawn from figure 5, one for each neuron
% figure 5 (with CSgauss stimulus only) - across-neurons distribution of cue peak values
