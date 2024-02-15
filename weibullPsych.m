function p = weibullPsych(S,valuesToTest,epsilon,steepnes)

p = 0.5+(1-0.5-epsilon/2).*(1-exp(-(valuesToTest./(1./10.^S)).^steepnes));
