function [prior] = priorCSF(PRIOR_GUESS,PARAMS_RANGE,plotPriors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [prior] = priorCSF(PRIOR_GUESS,PARAMS_RANGE)
% Definition of 4-D prior for the truncated log-parabola model of the csf.
% The marginal priors (for each parameter) are mostly flat uninformaitve
% defined by hyperbollic secant function, log-symmetric around mean.mode
% guees values defined in the PRIOR_GUESS structure:
% PRIOR_GUESS.f_max        peak spatial frequency, in cicles per degree
% PRIOR_GUESS.gamma_max    peak sensitivity (sensitivity = 1/contrast)
% PRIOR_GUESS.beta         bandwidthm that is the FWHM of the parabola (still not clear why this is defined as in 'octaves' unit according to Lesmes et al)
% PRIOR_GUESS.delta        truncation of the parabola in the log frequency range (in decimal log units)
% PRIOR_GUESS.confidence   
% the PARAMS.RANGE structure gives the range of the four parameters, which
% should be given uniformily spaced in log space
% It is important that the prior are symmetric in log space, otherwise the
% expected value of the posterior parameter values will not correspond with
% the mode, higuest values of the distribution
% Based on Lesmes et al. (2010), J Vis, 10(3):17 and Lu et al. (2023), Sci
% rep 13:16795
% JPO, Hamburg, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


priorFreq  = sech(PRIOR_GUESS.confidence.*(log10(PARAMS_RANGE.freqRange)-log10(PRIOR_GUESS.f_max)));
priorGamma = sech(PRIOR_GUESS.confidence.*(log10(PARAMS_RANGE.gammaRange)-log10(PRIOR_GUESS.gamma_max)));
priordelta = sech(PRIOR_GUESS.confidence.*(log10(PARAMS_RANGE.deltaRange)-log10(PRIOR_GUESS.delta)));
priorbeta  = sech(PRIOR_GUESS.confidence.*(log10(PARAMS_RANGE.betaRange)-log10(PRIOR_GUESS.beta)));

prior   = repmat(priorFreq',[1 PARAMS_RANGE.gammaN PARAMS_RANGE.deltaN PARAMS_RANGE.betaN]).*...
             repmat(priorGamma,[PARAMS_RANGE.freqN 1 PARAMS_RANGE.deltaN PARAMS_RANGE.betaN]).*...
             repmat(reshape(priordelta,[1 1 length(priordelta),1]),[PARAMS_RANGE.freqN PARAMS_RANGE.gammaN 1 PARAMS_RANGE.betaN]).*...
             repmat(reshape(priorbeta,[1 1 1 length(priorbeta)]),[PARAMS_RANGE.freqN PARAMS_RANGE.gammaN PARAMS_RANGE.deltaN 1]);

prior  = prior./sum(prior(:)); 

if plotPriors
    figure
    subplot(2,2,1)
    plot(log10(PARAMS_RANGE.freqRange) ,squeeze(sum(sum(sum(prior,2),3),4)))
    xlabel('peak frequency')
    subplot(2,2,2)
    plot(log10(PARAMS_RANGE.gammaRange),squeeze(sum(sum(sum(prior,1),3),4)))
    xlabel('gamma/ peak gain')
    subplot(2,2,3)
    plot(log10(PARAMS_RANGE.deltaRange) ,squeeze(sum(sum(sum(prior,1),2),4)))
    xlabel('delta, truncation')
    subplot(2,2,4)
    plot(log10(PARAMS_RANGE.betaRange) ,squeeze(sum(sum(sum(prior,1),2),3)))
    xlabel('beya, bandwidth')
end