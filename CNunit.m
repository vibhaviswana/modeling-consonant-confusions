function outp = CNunit(PSTHmultichan,list_CFs,indCF,psthbinwidth)

% Simulates a simple cochlear nucleus (CN) unit at a particular characteristic frequency (CF)
% Copyright 2021 Vibha Viswanathan. All rights reserved. 

% INPUTS:
% PSTHmultichan: peri-stimulus time histograms (PSTHs), size: CFs x PSTH bins
% list_CFs: list of CFs, size: 1 x CFs  
% indCF: index of the CF of this particular CN unit in list_CFs, size: scalar
% psthbinwidth: PSTH bin width (note: sampling rate = 1/psthbinwidth), size: scalar

% OUTPUT:
% outp: output PSTH of the CN unit, size: 1 x (PSTH bins + 1)

% select inhibitory receptive field (two octaves below & one octave above CF)
CF = list_CFs(indCF);
indAbove = ((list_CFs >= CF) & (list_CFs < CF*2)); % note: the inhibition is wideband, not lateral
indBelow = ((list_CFs < CF) & (list_CFs > CF/4));
indInhib = indAbove | indBelow;
numInhibConns = sum(indInhib); % no. inhib. connections
PSTHinhib = PSTHmultichan(indInhib,:);

% define time vector
t = 0:psthbinwidth:(size(PSTHmultichan,2)*psthbinwidth-psthbinwidth);

% excitatory synapse imposes low-pass filtering 
a = 1;
tau1 = 5e-3;
h_excit = exp(-t/tau1); % impulse response of excitatory synapse
h_excit = a*h_excit/sum(abs(h_excit));

% inhibitory synapse imposes low-pass filtering, delay, and inversion
b = -a/1.75;
delayInh = 2e-3; 
tau2 = 1e-3; 
h_inhib = exp(-(t-delayInh)/tau2).*(1+sign(t-delayInh))/2; % impulse resp. of inhibitory synapse
h_inhib = b*h_inhib/sum(abs(h_inhib));

% compute weighted sum of PSTH (firing rate/probability of firing) 
% output from excitatory and inhibitory synapses
c = 1/numInhibConns;
temp = max(conv(c*h_excit,sum(PSTHinhib,1),'full'),0);
outp_inhib = conv(h_inhib,temp,'full');
outp_excit = conv(h_excit,squeeze(PSTHmultichan(indCF,:)),'full');
outp = outp_inhib(1:numel(outp_excit)) + outp_excit;
outp = max(outp,0); % half-wave rectify
outp = outp(1:(size(PSTHmultichan,2)+1));

end

