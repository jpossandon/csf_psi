function p = weibullPsych(S,valuesToTest,epsilon,steepnes,chance)

p = chance+(1-chance-epsilon/2).*(1-exp(-(valuesToTest./(1./10.^S)).^steepnes));
