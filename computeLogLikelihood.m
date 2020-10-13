function llf=computeLogLikelihood(a,y,W,signal)
% function llf=computeLogLikelihood(a,y,signal)
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    Gamma=sqrt(signal.P)*W*a;   
    llf=-norm(y-(Gamma'*y)/(norm(Gamma)^2)*Gamma)^2/signal.sigma2; 