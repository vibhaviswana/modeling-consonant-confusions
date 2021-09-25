% Run auditory-nerve (AN) experiment using the Bruce et al. (2018) model  

% Copyright 2021 Vibha Viswanathan. All rights reserved. 

%% Setup

setup;

%% Run Bruce et al. (2018) model

psth = zeros(nconds,ntalkers,nconson,nCFs,numPSTHbins);
% SiDCmod is dichotic, so we need to store separate responses for the left and 
% right channels
SiSNM_psth_l = zeros(ntalkers,nconson,nCFs,numPSTHbins); 
SiSNM_psth_r = zeros(ntalkers,nconson,nCFs,numPSTHbins);
for m = 1:nconds
    for k = 1:ntalkers
        for p = 1:nconson
            disp(strcat('cond = ',num2str(m),', talker = ',num2str(k),', conson = ',num2str(p)));
            condi = conditions{m};
            cons = consonants{p};
            talk = talkers{k};
            [stim,fs] = audioread(fullfile(stimPath,strcat('stim_',condi,'_',...
                cons,'_Talker',talk,'.wav')));
            stim = resample(stim, fs_model/1000, round(fs_STeVI/1000)); % Resample to fs_model
            % Rescale stim to Pa units so that target has targetSPL. 20*log10(x/20e-6)=targetSPL
            stim = db2mag(targetSPL)*(20e-6)*stim/targRMS;
            % We have now converted all our stimuli to Pascal.
            if ((m~=1) && (m~=6)) % neither SiQuiet nor vocoded SiQuiet
                % speech-in-noise stimuli were created with speech starting 1 sec after masker
                stim = stim((1*fs_model+1):end,:);
            else % either SiQuiet or vocoded SiQuiet
                stim = stim;
            end
            stim = stim(1:T*fs_model,:); 
            inputl = stim(:,1)'; % select left stimulus channel
            if (strcmp(condi,'SiSNM_-18dB')) % process both channels
                inputr = stim(:,2)'; % select right stimulus channel
            end
            % Run Bruce et al. (2018) model
            for whichCF = 1:nCFs
                CF = list_CFs(whichCF);
                % left audio
                vihc = model_IHC_BEZ2018(inputl,CF,nrep,1/fs_model,...
                    T,cohc,cihc,species);
                [out1,~,~] = model_Synapse_BEZ2018(vihc,...
                    CF,nrep,1/fs_model,noiseType,implnt,spont,tabs,trel);
                temp_psthl = sum(reshape(out1,...
                    nTimePtsPerPSTHbin,numPSTHbins)); % Calculate proper PSTH (unit: spikes/s)
                temp_psthl = temp_psthl/nrep/psthbinwidth;
                psth(m,k,p,whichCF,:) = temp_psthl;
                if (strcmp(condi,'SiSNM_-18dB')) % process both channels
                    % right audio
                    vihc = model_IHC_BEZ2018(inputr,CF,nrep,1/fs_model,...
                        T,cohc,cihc,species);
                    [out1,~,~] = model_Synapse_BEZ2018(vihc,...
                        CF,nrep,1/fs_model,noiseType,implnt,spont,tabs,trel);
                    temp_psthr = sum(reshape(out1,...
                        nTimePtsPerPSTHbin,numPSTHbins)); % Calculate proper PSTH (unit: spikes/s)
                    temp_psthr = temp_psthr/nrep/psthbinwidth;
                    psth(m,k,p,whichCF,:) = temp_psthl + temp_psthr;
                    SiSNM_psth_l(k,p,whichCF,:) = temp_psthl;
                    SiSNM_psth_r(k,p,whichCF,:) = temp_psthr;
                end
            end
        end
    end
end

save('Bruce2018ANoutput_CVstim.mat','psth','SiSNM_psth_l','SiSNM_psth_r');

