function [nextFreqToTest,nextContrastToTest] = csf_stim_select(priorCSF,freqToTest,contrastToTest,PARAMS_RANGE,PARAMS_SELECT)

NSamples = PARAMS_SELECT.NSamples ;

% select a random sample size n from prior distribution parameters
x                 = discretesample(priorCSF(:), NSamples);
% get the respective indexes for fmax, gamma, delta and beta
[Fix,Gix,Dix,Bix] = ind2sub(size(priorCSF),x);

for ff = 1:length(freqToTest)
    for cc = 1:length(contrastToTest)
        
        % calcualte weibull psychometric function for given stimulus for all
        % sampled S
        for nn = 1:NSamples
            S = csf(PARAMS_RANGE.freqRange(Fix(nn)),...
                PARAMS_RANGE.gammaRange(Gix(nn)),...
                PARAMS_RANGE.deltaRange(Dix(nn)),...
                PARAMS_RANGE.betaRange(Bix(nn)),...
                freqToTest(ff));
            p(nn) = weibullPsych(S,contrastToTest(cc),PARAMS_SELECT.epsilon,PARAMS_SELECT.steepnes,PARAMS_SELECT.chance);
        end

        % Calculate information gain
            avgEntropies        =  sum(-p.*log(p)-((1-p).*log(1-p)))./NSamples;
            avgPCorrect         = sum(p)./NSamples;
            entropy_avgPCorrect = sum(-avgPCorrect.*log(avgPCorrect)-((1-avgPCorrect).*log(1-avgPCorrect)));
            expInfoGain(ff,cc)  = entropy_avgPCorrect - avgEntropies;
    end
end

% select one stimuli from top quantile (according to parameter unifDecile)
[~,eIGsort]         = sort(expInfoGain(:),'descend');
SampleTopQuantile   = randperm(round(numel(expInfoGain).*PARAMS_SELECT.unifDecile),1);
[nextFreqToTest,nextContrastToTest] = ind2sub([length(freqToTest) length(contrastToTest)],eIGsort(SampleTopQuantile));

nextFreqToTest      = freqToTest(nextFreqToTest);
nextContrastToTest  = contrastToTest(nextContrastToTest);