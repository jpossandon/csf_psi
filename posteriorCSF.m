function posterior = posteriorCSF(prior,freqTested,contrastTested,correct,PARAMS_RANGE,PARAMS_SELECT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% posterior = posteriorCSF(PARAMS_RANGE,freqTested,contrastTested,correct)
% calculates posterior according to Bayes Rule and a trial answer
% jpo, 2024, Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all possible combination of parameters that are present in the prior
try % combinations was introduced for Matlab R2023a
    T = combinations( PARAMS_RANGE.freqRange, PARAMS_RANGE.gammaRange,PARAMS_RANGE.deltaRange, PARAMS_RANGE.betaRange);
    T.Properties.VariableNames = {'p_f','gamma','delta','bw'};
catch
    [an, bn, cn, dn] = ndgrid(PARAMS_RANGE.freqRange, PARAMS_RANGE.gammaRange,PARAMS_RANGE.deltaRange, PARAMS_RANGE.betaRange); 
    T = table(an(:), bn(:), cn(:), dn(:),'VariableNames',{'p_f','gamma','delta','bw'});
    clear an bn cn dn
end

% calculate the probability of corect or incorrect response for all the CSF
% parameter space and the given frequency and contrast tested 
S       = csf(T.p_f,T.gamma,T.delta,T.bw,freqTested);
pOBS    = weibullPsych(S,contrastTested,PARAMS_SELECT.epsilon,PARAMS_SELECT.steepnes);

if correct == 0
    pOBS = 1-pOBS;
end

% reshape pOBS to fit the prior
[lia1,loc_p_f]      = ismember(T.p_f,PARAMS_RANGE.freqRange);
[lia1,loc_gamma]    = ismember(T.gamma,PARAMS_RANGE.gammaRange);
[lia1,loc_delta]    = ismember(T.delta,PARAMS_RANGE.deltaRange);
[lia1,loc_beta]     = ismember(T.bw,PARAMS_RANGE.betaRange);

INDEX               = sub2ind(size(prior),loc_p_f,loc_gamma,loc_delta,loc_beta);
clear lia1 loc_p_f loc_gamma loc_beta
pOBS2               = nan(size(prior));
pOBS2(INDEX)        = pOBS;

% calculate the total probability of the data (correct or incorrect) across
% the complete prior space
pData               = sum(pOBS2(:).*prior(:));
% calcualte the posterior,
posterior           = prior.*pOBS2./pData;






% % slower loop alternative to test (for a 40x60x40x27 prior the code above takes 200 ms and the loop one 23s)
% tic
% pOBStest = nan(size(prior));
% for f=1:length(PARAMS_RANGE.freqRange)
%     for g = 1:length(PARAMS_RANGE.gammaRange)
%         for d = 1:length(PARAMS_RANGE.deltaRange)
%             for b = 1:length(PARAMS_RANGE.betaRange)
%                  S = csf(PARAMS_RANGE.freqRange(f),...
%                 PARAMS_RANGE.gammaRange(g),...
%                 PARAMS_RANGE.deltaRange(d),...
%                 PARAMS_RANGE.betaRange(b),...
%                 nextFreqToTest);
%                 pOBStest(f,g,d,b) = weibullPsych(S,nextContrastToTest,PARAMS_SELECT.epsilon,PARAMS_SELECT.steepnes);
%             end
%         end
%     end
% end
% toc