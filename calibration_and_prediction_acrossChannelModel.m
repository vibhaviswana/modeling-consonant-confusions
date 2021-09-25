% Calibrates the across-channel model with SiSSN data, 
% then predicts confusions in novel conditions.
% Copyright 2021 Vibha Viswanathan. All rights reserved.  

setup;

%% Process model output (neural confusion matrices)

temp = load(fullfile(ANsimResultPath,'acrossChannelModelOutput'));
neuralConf = temp.confMat; 
neuralConf_temp = zeros(nconds, ntalkers, nconsonants, nconsonants);
wts = ones(nCFs, nModBands); 
for c = 1:nconds
    for t = 1:ntalkers
        for k = 1:nconsonants
            for l = 1:nconsonants
                temp = squeeze(neuralConf(c, t, :, :, k, l));
                temp = temp .* wts;
                denr = numel(temp);
                numr = ceil(1.0 * denr); 
                temptop = maxk(temp(:), numr);
                neuralConf_temp(c, t, k, l) = mean(temptop);
            end
        end
    end
end

neuralConf_temp = squeeze(mean(neuralConf_temp, 2));
neuralConf = neuralConf_temp;
neuralConf = neuralConf./...
    repmat(sum(neuralConf,3),...
    1,1,nconsonants); % divide each row by sum of all elements so that you get a proportion

% Extract SiSSN neural confusions
SiSSN_neuralConf = squeeze(neuralConf(2,:,:));

%% Process measured perceptual confusion matrices

fname1 = 'preprocessed_CI_RepExpt1&2_Test_results.csv'; 
% Read results as dataframe with 6 columns
% Column 1: Subject ID
% Column 2: Talker
% Column 3: Condition
% Column 4: consonant presented (ranges from 1-20)
% Column 5: consonant heard/reported (ranges from 1-20)
% Column 6: Whether or not subjects answered correctly
% Column 7: Response time
% Column 8: Position in confusion matrix
dTable1 = readtable(fullfile(perceptualDatapath,fname1));
uniqueSubjTalkerCombos = unique(dTable1(:,[1,2]));
nustc1 = size(uniqueSubjTalkerCombos,1);

confusionMatrix1 = zeros(nconds,nconsonants,nconsonants);
for c = 1:nconds
    indc = find(dTable1.Condition == c);
    Tcond = dTable1(indc,:);
    for i = 1:nconsonants
        index3 = find(Tcond.ConsPresented == i);
        responses = Tcond.ConsHeard(index3);
        for j = 1:numel(responses)
            confusionMatrix1(c,i,responses(j)) = ...
                confusionMatrix1(c,i,responses(j)) + 1;
        end
    end
end

% DO NOT normalize confusion matrices to a fixed overall intelligibility as
% we want to predict intelligibility too.
% Instead normalize each confusion matrix by the max value any element can take
confusionMatrix1 = confusionMatrix1/nustc1;

SiSSN_perceptConf = squeeze(confusionMatrix1(2,:,:)); % Perceptual SiSSN confusions

%% Calibration with SiSSN
% Use SiSSN data for calibration. Separately calibrate each model.

% Similarity between neural and perceptual metrics for SiSSN
corr_SiSSN = corrcoef(SiSSN_neuralConf(:),SiSSN_perceptConf(:));

% Fit a sigmoid/logistic function mapping neural to perceptual confusions
% using SiSSN
sigmoid = @(a, b, x) 1. / (1 + exp(-b*(x-a))); % y = 1/(1 + e ^ (-k(x-x0)))

fitobject = fit(SiSSN_neuralConf(:),SiSSN_perceptConf(:),sigmoid);
ylogis = sigmoid(fitobject.a, fitobject.b, SiSSN_neuralConf(:));
figure;
[~,inds] = sort(SiSSN_neuralConf(:),'ascend');
plot(SiSSN_neuralConf(inds),ylogis(inds),'-','linewidth',2); hold on;
plot(SiSSN_neuralConf(:),SiSSN_perceptConf(:),'o','markersize',8);
set(gca,'FontSize',16);
legend({'Fit','Data'},'location','southeast','fontsize',18);
xlabel('Neural confusion matrix entry for across-channel model','fontsize',22);
ylabel('Perceptual confusion matrix entry','fontsize',22);

% Predict confusions for other conditions
predConds = 3; % choose one condition to predict
predictionAccuracy = zeros(1,numel(predConds));
for m = 1:numel(predConds)
    k = predConds(m);
    neuralconfmat = squeeze(neuralConf(k,:,:));
    oinds = true(size(neuralconfmat)); % use all entries
    neuralconfmat = neuralconfmat(oinds); 
    predicted_perceptconfmat = sigmoid(fitobject.a, fitobject.b, neuralconfmat(:));
    perceptconfmat = squeeze(confusionMatrix1(k,:,:));
    perceptconfmat = perceptconfmat(oinds);
    temp_corr = corrcoef(perceptconfmat(:),predicted_perceptconfmat(:));
    predictionAccuracy(m) = temp_corr(1,2);
end
 
%% Save predicted and actual for the above-selected condition
A = reshape(predicted_perceptconfmat, nconsonants, nconsonants);  % predicted confusions
B = reshape(perceptconfmat, nconsonants, nconsonants); % measured confusions
save predictions_SiB_acrossChannel A B; % Change filename to match predConds


