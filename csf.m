function [S] = csf(p_f,gamma,delta,bw,spatFreq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [S] = csf(gamma,p_f,bw,delta,spatFreq)
% Generates log10(sensitivity) values at spatFreqs, according to the truncated log
% parabola model, defined by four parameters
% gamma  - peak gain/sensitivity; that is 1/contrast, with contrast in values between
%           0 and 1 
% p_f    - peak spatial frequency in cycles per degree (where the parabola peaks)
% bw     - bandwidth of the parabola, in 'octaves
% delta  - truncation level of the low spatial frequency side, in decimal
%           log units below the peak
%
% Based on Lesmes et al. (2010), J Vis, 10(3):17 and Lu et al. (2023), Sci
% rep 13:16795
% JPO, 2024, Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

kappa   = 4./log10(2);
S    = log10(gamma)-kappa.*((log10(spatFreq./p_f))./bw).^2; %log parabola

if length(spatFreq)>=1 & length(p_f)==1
    S(spatFreq <p_f & S<log10(gamma)-delta) = log10(gamma)-delta; % truncation
elseif length(spatFreq)==1 & length(p_f)>1
    spatFreq = repmat(spatFreq,length(p_f),1);
    whichToTtruncate = spatFreq < p_f & S<log10(gamma)-delta;
    S(whichToTtruncate) = log10(gamma(whichToTtruncate))-delta(whichToTtruncate); % truncation
end