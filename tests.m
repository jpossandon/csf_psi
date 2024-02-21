%%
% Contrast is here defined betwenn 0 and 1. For periodic stimuli like gratings/sinusoid the
% range of 0 to 1 corresponds with Michelson's contrast for stimuli that
% vary around the middle gray value the screen show ('127'). In a 8-bit monitor, and considering a
% luminnance value of 0 as no light and without considering gamma correction
% 'luminance' values can go between 1 ((255-0)/(255+0)) and 0.0079
% (128-126)/(128+126). (See Tyler, 1992, bit-stealing techinche for a way
% to increase the possible contrast values in an 8bit monitor)

% Frequency is defined as cycles per degrees, and thus determined by the
% monitor resolution, size and the distnace of the observer

% SET PARAMETERS

% STIMULUS SET
% It does not necessarily need to be as the CSF prior parameter, both frequencies and contrasts 
% to test have to be possible in the setup

freqToTest      = logspace(log10(.2),log10(36),20);

contrastToTest  = logspace(log10(.002),log10(1),30); % this should be adjusted 

% CSF PRIOR DEFINITION
% The CSF is definned by a log-parabola model with truncation as explained
% by Lesmes et al. (2010). See csf.m for more detailed definitions

% Guesses, this should be adjusted for special population. For
% sight-recovery participant follow Kalia et al, 2014
PARAMS.PRIOR_GUESS.f_max        = 2.5 ;  % peak spatial frequency, in cicles per degree
PARAMS.PRIOR_GUESS.gamma_max    = 100 ;  % peak sensitivity (sensitivity = 1/contrast)
PARAMS.PRIOR_GUESS.beta         = 3 ;    % bandwidthm that is the FWHM of the parabola (still not clear why this is defined as in 'octaves' unit according to Lesmes et al)
PARAMS.PRIOR_GUESS.delta        = 0.5 ;  % truncation of the parabola in the log frequency range (in decimal log units)
PARAMS.PRIOR_GUESS.confidence   = 1;

% range of the prior parameters (PRIORS NEED TO BE SYMMETRIC!!)
PARAMS.RANGE.gammaN         = 40;
PARAMS.RANGE.deltaN         = 40;
PARAMS.RANGE.betaN          = 27;
PARAMS.RANGE.freqN          = 40;
PARAMS.RANGE.freqRange      = logspace(log10(.2),log10(36),PARAMS.RANGE.freqN);
PARAMS.RANGE.gammaRange     = logspace(log10(10),log10(1000),PARAMS.RANGE.gammaN); 
PARAMS.RANGE.deltaRange     = logspace(log10(0.05),log10(5),PARAMS.RANGE.deltaN);
PARAMS.RANGE.betaRange      = logspace(log10(1),log10(9),PARAMS.RANGE.betaN);

% how to select next stimuli parameters
PARAMS.SELECT.epsilon     = 0.04;   % stimulus independe error or lapses 1-epsilon is the max performance
PARAMS.SELECT.chance      = 0.25;   % lower bound of the psychometric function, in a n-AFC task correspond to 1/n
PARAMS.SELECT.steepnes    = 3;      % the steepnes/slope of the Weibull psychometric function, it is assumed to be similar for all frequencies and observers, which might be wrong, values used are between 2 and 3.5
PARAMS.SELECT.NSamples    = 500;    % sampling of the posterior to evaluate next stimuli
PARAMS.SELECT.unifDecile  = .1;     % instead of selecting the stimuli that gives the max expected information gain, the slection is a random selection between the unifDecile of stimuli (this is just mentioned in Lesmes et al, p14)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment simulation

% starting parameters prior
plotPriors  = 1;    % check that prior are symmetric
[prior]     = priorCSF(PARAMS.PRIOR_GUESS,PARAMS.RANGE,plotPriors);

% the csf to find in this 'experiment'
gamma       = 100;
p_f         = 2;
bw          = 1.5;
delta       = 1;
[S]         = csf(p_f,gamma,delta,bw,PARAMS.RANGE.freqRange);

% figure to check the simulation
figure
% plot of the csf to find (in black) and marginal 2-D prior mass functions
% In each iteration a red/black dot appear at the frequency and contrast tested, a csf 
% according to the expected value of the parameters given the posterior is plotted in light
% green, and the posterior mass functions are updated. The last csf (i.e., the result of the experiment)
% is shown in red at the end of the simulation

subplot('Position',[.1 .25 .5 .5])
ax1 = plotCSF(PARAMS.RANGE.freqRange,S); hold on
title(sprintf('pf:%1.1f cpd  gamma:%1.1f delta:%1.1f bw:%1.1f',p_f,gamma,delta,bw))
ttext = text(log2(18),log10(100),sprintf('trial '));
subplot('Position',[.7 .6 .25 .25])
p_f_vs_gamma = squeeze(sum(sum(prior,3),4));
% pcolor(log10(PARAMS.RANGE.gammaRange),log2(PARAMS.RANGE.freqRange),p_f_vs_gamma)
imagesc(log10(PARAMS.RANGE.gammaRange),log2(PARAMS.RANGE.freqRange),p_f_vs_gamma),hold on
xlabel('Peak Sensitivity')
ylabel('Peak Freq.')
title('posterior density')
ax2 = gca;
yticks2 = [.5 1 2 4 8 16 32];
xticks2 = [1 5 10 20 40 80 160 320];
set(ax2,'XTick',log10(xticks2),'XTickLabel',xticks2,...
    'YTick',log2(yticks2),'YTickLabel',yticks2)

subplot('Position',[.7 .1 .25 .25])
delta_vs_bw = squeeze(sum(sum(prior,1),2));

imagesc(log10(PARAMS.RANGE.betaRange),log10(PARAMS.RANGE.deltaRange),delta_vs_bw),hold on
ylabel('Truncation')
xlabel('Bandwidth')
title('posterior density')
ax3 = gca;
xticks3 = [1 2 4 8];
yticks3 = [.1 .5 1 2 4];
set(ax3,'XTick',log10(xticks3),'XTickLabel',xticks3,...
    'YTick',log10(yticks3),'YTickLabel',yticks3)

%%
nTrials = 70; % number of trials
lcolor = [0 1 0 .1];

for tt = 1:nTrials

    % selection new stimuli to test
    [nextFreqToTest,nextContrastToTest] = csf_stim_select(prior,freqToTest,contrastToTest,PARAMS.RANGE,PARAMS.SELECT);
   
    % correct or incorrect according to the csf defined above and the
    % weibull psychometric function
    [St]        = csf(p_f,gamma,delta,bw,nextFreqToTest);
    p           = weibullPsych(St,nextContrastToTest,PARAMS.SELECT.epsilon ,PARAMS.SELECT.steepnes,PARAMS.SELECT.chance);
    correct     = rand(1)<p;
     
    % update prior
    prior       = posteriorCSF(prior,nextFreqToTest,nextContrastToTest,correct,PARAMS.RANGE,PARAMS.SELECT);

    % this iteration estimate of each parameter value (the expected value of the prior marginals)
    p_f_est     = 10.^(sum(sum(sum(sum(prior,2),3),4).*log10(PARAMS.RANGE.freqRange')));
    gamma_est   = 10.^sum(sum(sum(sum(prior,1),3),4).*log10(PARAMS.RANGE.gammaRange));
    delta_est   = 10.^sum(squeeze(sum(sum(sum(prior,1),2),4)).*log10(PARAMS.RANGE.deltaRange)');
    beta_est    = 10.^sum(squeeze(sum(sum(sum(prior,1),2),3)).*log10(PARAMS.RANGE.betaRange)');
   
    % plot of the csf according to the estimates and the changes in the
    % posterior
    [Stt] = csf(p_f_est,gamma_est,delta_est,beta_est,PARAMS.RANGE.freqRange);
    axes(ax1)
    if tt == nTrials, lcolor = [1 0 0];end
    plot(log2(PARAMS.RANGE.freqRange),Stt,'-','Color',lcolor,'LineWidth',1); 
    ttext.String = sprintf('trial %d',tt);
    scatter(log2(nextFreqToTest),log10(1./nextContrastToTest),17,'o',...
        'MarkerEdgeColor',[correct 0 0],'MarkerFaceColor',[correct 0 0],'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
    axes(ax2)
    p_f_vs_gamma = squeeze(sum(sum(prior,3),4));
    imagesc(log10(PARAMS.RANGE.gammaRange),log2(PARAMS.RANGE.freqRange),p_f_vs_gamma)

    axes(ax3)
    delta_vs_bw = squeeze(sum(sum(prior,1),2));
    imagesc(log10(PARAMS.RANGE.betaRange),log10(PARAMS.RANGE.deltaRange),delta_vs_bw)
end

