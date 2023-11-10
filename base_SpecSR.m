function [HHmean,HHstd] = base_SpecSR(StimF,StimA,bipolar,dt_a,pw,SR,Trials,Time)

%% Outputs
% HHmean and HHstd represent the mean and standard deviation of the
% half harmonic/stimulation frequency aliasing artefact over all trials

%Properties of Pulse train
%% StimF is the stimulation frequency of the pulse train
if isempty(StimF)
    StimF = 130; % 130 chosen as what was for figure in publication and common DBS frequency
end
%% StimA is the stimulation amplitude of the pulse train, does not have much more effect than scaling PSD
if isempty(StimA)
    StimA = 10; % 10 was used for figure in publication
end
%% bipolar determines whether pulse train has single positive pulse (false) or bipolar pulses (true)
if isempty(bipolar)
    bipolar = false; % false was used for figure in publication
end
%% dt_a represents the time steps for the generated pulse train to represent an analogue signal.
%% Smaller is better
if isempty(dt_a)
    dt_a = 1*10^(-6); % 1*10^(-6) chosen as far less than sampling rates and no clear differences
%     seen at smaller time steps in PSD, 1*10^(-6) used in publication
end
%% pw is the pulse width of the pulse train, if bipolar then both pulses take same width of pw
if isempty(pw)
    pw = 90*10^(-6); % 90*10^(-6) was used for figure in publication and 90mus is common pulse width for DBS
end

%Property important for the ADC
%% SRVect represent the range of sampling rates chosen for figure plotting
if isempty(SR)
    SR = 1000; % 250:1:5000 was used for figure in publication and provides a wide range
end

%Properties of Plotting and Robustness
%% Trials is the simulations with different starting points that the PSD is calculated for.
%% This improves accuracy of the mean of simulations with the more Trials performed
if isempty(Trials)
    Trials = 100; % 10 was used for figure in publication
end
%% Time represents the length of a single simulation
%% The greater Time is the higher the overall accuracy will be for a single simulation
if isempty(Time)
    Time = 1; % 1 was used for figure in publication and captures enough pulse cycles for 130Hz
end

freqs = 0:0.1:2*StimF;
HalfHarm_rec = NaN(1,Trials);
Harm1_rec = NaN(1,Trials);

if bipolar
    Pulse = StimA.*[ones(round(pw/dt_a),1);-1.*ones(round(pw/dt_a),1)];
else
    Pulse = StimA.*ones(round(pw/dt_a),1);
end

[~,HH_ind] = min(abs(freqs - StimF/2));
[~,H1_ind] = min(abs(freqs - StimF));

SampT = 1/SR; Sampdt = SampT/dt_a;
for i_Tr = 1:Trials

    Ts = 1/StimF;
    RestP = zeros(round(Ts/dt_a)-length(Pulse),1);
    StimIt = [Pulse;RestP];
    StimSeq_a = repmat(StimIt,ceil(Time/Ts),1);

    iniSamp = round(((Sampdt-1)*rand(1))+1);
    SampPoints = round(iniSamp+((0:(Time*SR)-1).*Sampdt));
    SampPoints(SampPoints > length(StimSeq_a)) = [];
    StimSeq_d = StimSeq_a(SampPoints);

    PSD = pwelch(StimSeq_d-mean(StimSeq_d),length(StimSeq_d),[],freqs,SR);
    HalfHarm_rec(1,i_Tr) = PSD(HH_ind(1));
    Harm1_rec(1,i_Tr) = PSD(H1_ind(1));
end
Harm1_rec(Harm1_rec == 0) = 0.00000001;
NormHarm = HalfHarm_rec./Harm1_rec;

HHmean = mean(NormHarm);
HHstd = std(NormHarm);

if HHmean >= 0.1
    disp(['Half harmonic artefact is prominent = ' num2str(HHmean) ', choose new sampling rate'])
elseif HHmean + HHstd >= 0.1
    disp(['Half harmonic artefact is frequently prominent = ' num2str(HHmean) ', choose new sampling rate'])
elseif HHmean + HHstd >= 0.01
    disp(['Half harmonic artefact is low but present = ' num2str(HHmean) ', be careful with this sampling rate'])
else
    disp(['Half harmonic artefact should not impact recording = ' num2str(HHmean)])
end