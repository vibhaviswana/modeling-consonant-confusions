% Across-channel temporal-coherence processing model (based on computations
% that exist in the cochlear nucleus) that accounts for across-channel
% modulation interference to predict consonant confusions.

% Copyright 2021 Vibha Viswanathan. All rights reserved.

%% Setup

setup;

% Load PSTHs from AN model experiment
f1 = load(fullfile(ANsimResultPath,'Bruce2018ANoutput_CVstim.mat'));
neuralResp = f1.psth;
SiSNM_respl = f1.SiSNM_psth_l;
SiSNM_respr = f1.SiSNM_psth_r;

%% Pass nerve outputs to cochlear nucleus (CN) unit

Cdur_bins = Cdur_bins + 50; % added for across-channel model only

cn_resp = zeros(nconds,ntalkers,nconson,nCFs,Cdur_bins+1);
cn_SiSNM_l = zeros(ntalkers,nconson,nCFs,Cdur_bins+1);
cn_SiSNM_r = zeros(ntalkers,nconson,nCFs,Cdur_bins+1);

for aTalker = 1:ntalkers
    for aCons = 1:nconson
        consBins = consBinsStart(aTalker,aCons):(consBinsStart(aTalker,aCons)+Cdur_bins);
        for aCond = 1:nconds
            resp = squeeze(neuralResp(aCond,aTalker,aCons,:,consBins)); % CF x time matrix
            for p = 1:nCFs
                CF_hertz = list_CFs(p);
                cn_resp(aCond,aTalker,aCons,p,:) = CNunit(resp,list_CFs,p,psthbinwidth);
            end
            if (aCond == 4) % SiDCmod
                respl = squeeze(SiSNM_respl(aTalker,aCons,:,consBins)); % CF x time matrix
                respr = squeeze(SiSNM_respr(aTalker,aCons,:,consBins)); % CF x time matrix
                for p = 1:nCFs
                    CF_hertz = list_CFs(p);
                    cn_SiSNM_l(aTalker,aCons,p,:) = CNunit(respl,list_CFs,p,psthbinwidth);
                    cn_SiSNM_r(aTalker,aCons,p,:) = CNunit(respr,list_CFs,p,psthbinwidth);
                end
            end
        end
    end
end

%% Dynamic time warping 

consPairs = zeros(2,nconson*nconson);
k = 0;
for i = 1:nconson
    for j = 1:nconson
        k = k + 1;
        consPairs(1,k) = i;
        consPairs(2,k) = j;
    end
end

% Compute dynamic time warping for consonants in quiet
for aTalker = 1:ntalkers
    for aConsPair = 1:(nconson*nconson)
        aConsHeard = consPairs(1,aConsPair); % consonant reported
        aConsPres = consPairs(2,aConsPair); % consonant presented
        envQuiet = squeeze(cn_resp(1,aTalker,aConsHeard,:,:));
        envCond = squeeze(cn_resp(1,aTalker,aConsPres,:,:));
        for aCF = 1:nCFs
            ord = 0;
            framelen = 1;
            [~,tx,ty] = dtw(sgolayfilt(envQuiet(aCF,:), ord, framelen),...
                sgolayfilt(envCond(aCF,:), ord, framelen));
            ix{aTalker,aConsPair,aCF} = tx;
            iy{aTalker,aConsPair,aCF} = ty;
        end
    end
end

%% Modulation filterbank per Relano-Iborra et al. (2016)

env_byModCF = zeros(nconds,ntalkers,nCFs,nconson,nModBands,Cdur_bins+1-1);
env_byModCF_SiDCmod_l = zeros(ntalkers,nCFs,nconson,nModBands,Cdur_bins+1-1);
env_byModCF_SiDCmod_r = zeros(ntalkers,nCFs,nconson,nModBands,Cdur_bins+1-1);
for aCond = 1:nconds
    for aTalker = 1:ntalkers
        for aCF = 1:nCFs
            for aCons = 1:nconson
                inp = squeeze(cn_resp(aCond,aTalker,aCons,aCF,:));
                outp = modFbank_v3(inp,fs_PSTH,cf_mod);
                env_byModCF(aCond,aTalker,aCF,aCons,:,:) = outp;
                if (aCond == 4) % SiDCmod
                    inpl = squeeze(cn_SiSNM_l(aTalker,aCons,aCF,:));
                    outpl = modFbank_v3(inpl,fs_PSTH,cf_mod);
                    env_byModCF_SiDCmod_l(aTalker,aCF,aCons,:,:) = outpl;
                    inpr = squeeze(cn_SiSNM_r(aTalker,aCons,aCF,:));
                    outpr = modFbank_v3(inpr,fs_PSTH,cf_mod);
                    env_byModCF_SiDCmod_r(aTalker,aCF,aCons,:,:) = outpr;
                end
            end
        end
    end
end

%% Derive neural confusion matrices

confMat = zeros(nconds,ntalkers,nCFs,nModBands,nconson,nconson);

for aCond = 1:nconds
    for aTalker = 1:ntalkers
        for aConsPair = 1:(nconson*nconson)
            aConsHeard = consPairs(1,aConsPair); % consonant reported
            aConsPres = consPairs(2,aConsPair); % consonant presented
            for aModBand = 1:nModBands
                envQuiet = squeeze(env_byModCF(1,aTalker,:,aConsHeard,aModBand,:));
                if (aCond == 4) % SiDCmod
                    envCond_l = squeeze(env_byModCF_SiDCmod_l(aTalker,:,aConsPres,aModBand,:));
                    envCond_r = squeeze(env_byModCF_SiDCmod_r(aTalker,:,aConsPres,aModBand,:));
                    for aCF = 1:nCFs
                        ntimes = Cdur_bins+1-1;
                        ix_adj = ix{aTalker,aConsPair,aCF};
                        ix_adj(ix_adj > ntimes) = ntimes;
                        iy_adj = iy{aTalker,aConsPair,aCF};
                        iy_adj(iy_adj > ntimes) = ntimes;
                        warpedX = envQuiet(aCF, ix_adj);
                        warpedY_l = envCond_l(aCF, iy_adj);
                        warpedY_r = envCond_r(aCF, iy_adj);
                        corrinp1 = warpedX(:);
                        corrinp2_l = warpedY_l(:);
                        corrinp2_r = warpedY_r(:);
                        [temp_l,~] = corrcoef(corrinp1,corrinp2_l);
                        temp_l = temp_l(1,2);
                        inds = ((temp_l<0) | isnan(temp_l));
                        temp_l(inds) = 0;
                        temp_l = temp_l.^2; % convert to variance explained
                        [temp_r,~] = corrcoef(corrinp1,corrinp2_r);
                        temp_r = temp_r(1,2);
                        inds = ((temp_r<0) | isnan(temp_r));
                        temp_r(inds) = 0;
                        temp_r = temp_r.^2; % convert to variance explained
                        % Use "better" ear (i.e., use max instead of mean)
                        confMat(aCond,aTalker,aCF,aModBand,aConsPres,aConsHeard) = max(temp_l,temp_r);
                    end
                else % all other conditions presented diotic stimuli
                    envCond = squeeze(env_byModCF(aCond,aTalker,:,aConsPres,aModBand,:));
                    for aCF = 1:nCFs
                        ntimes = Cdur_bins+1-1;
                        ix_adj = ix{aTalker,aConsPair,aCF};
                        ix_adj(ix_adj > ntimes) = ntimes;
                        iy_adj = iy{aTalker,aConsPair,aCF};
                        iy_adj(iy_adj > ntimes) = ntimes;
                        warpedX = envQuiet(aCF, ix_adj);
                        warpedY = envCond(aCF, iy_adj);
                        corrinp1 = warpedX(:);
                        corrinp2 = warpedY(:);
                        [temp,~] = corrcoef(corrinp1,corrinp2);
                        temp = temp(1,2);
                        inds = ((temp<0) | isnan(temp));
                        temp(inds) = 0;
                        temp = temp.^2; % convert to variance explained
                        confMat(aCond,aTalker,aCF,aModBand,aConsPres,aConsHeard) = temp;
                    end
                end
            end
        end
    end
end

save('acrossChannelModelOutput','confMat');


