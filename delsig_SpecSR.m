function [HHmean,HHstd] = delsig_SpecSR(StimF,StimA,bipolar,dt_a,pw,Pfreq,Threshold,Trials,Time)

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
%% Pfreqs corresponding frequencies of predetermined period (PVect) that pulses are summed in (approx. sampling rate)
if isempty(Pfreq)
    Pfreq = 1000; % linspace(250,5000,10000) was used for figure in publication
end
%% Threshold the integral wave is required to cross to indicate recording of a pulse and resetting of intergral wave
if isempty(Threshold)
    Threshold = 900; % 900 was used for figure in publication and is large enough to trigger with StimA at 10
end

%Properties of Plotting and Robustness
%% Trials is the simulations with different starting points that the PSD is calculated for.
%% This improves accuracy of the mean of simulations with the more Trials performed
if isempty(Trials)
    Trials = 2; % 10 was used for figure in publication
end
%% Time represents the length of a single simulation
%% The greater Time is the higher the overall accuracy will be for a single simulation
if isempty(Time)
    Time = 0.1; % 1 was used for figure in publication and captures enough pulse cycles for 130Hz
end

freqs = 0:0.1:2*StimF;
PVect = round(1./Pfreq,6);

HalfHarm_rec = NaN(1,Trials);
Harm1_rec = NaN(1,Trials);

if bipolar
    Pulse = StimA.*[ones(round(pw/dt_a),1);-1.*ones(round(pw/dt_a),1)];
else
    Pulse = StimA.*ones(round(pw/dt_a),1);
end

Ts = 1/StimF;
RestP = zeros(round(Ts/dt_a)-length(Pulse),1);
StimIt = [Pulse;RestP];
StimSeq_a = repmat(StimIt,ceil(Time/Ts),1);

PulseInt = NaN(length(StimSeq_a),1);
PulseInt(1) = StimSeq_a(1);
for i_Pint = 2:length(StimSeq_a)
%     PulseInt(i_Pint) = sum(StimSeq_a(1:i_Pint));
    PulseInt(i_Pint) = PulseInt(i_Pint-1) + StimSeq_a(i_Pint);
end

[~,HH_ind] = min(abs(freqs - StimF/2));
[~,H1_ind] = min(abs(freqs - StimF));

P = PVect; Plength = round(P/dt_a); SR = 1/P;
for i_Tr = 1:Trials

    iniSamp = round((Plength-1)*rand(1))+1;
    NumP = floor((length(StimSeq_a)-iniSamp)/Plength);
    StSelection = PulseInt(iniSamp:(iniSamp + (NumP*Plength))) - PulseInt(iniSamp);
    Counter = NaN(1,NumP);
    for i_Counter = 1:NumP
        if ~bipolar
            Counter(i_Counter) = floor((StSelection(round(i_Counter*Plength)) -...
                StSelection(round(((i_Counter-1)*Plength)+1)))/Threshold);
            StSelection = StSelection - Counter(i_Counter)*Threshold;
        end
    end

%         PSD = pwelch(StimSeq_d-mean(StimSeq_d),length(StimSeq_d),[],freqs,SR);
    PSD = pwelch(Counter-mean(Counter),length(Counter),[],freqs,SR);
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