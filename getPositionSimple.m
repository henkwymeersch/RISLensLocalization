function [Pest,LLF]=getPositionSimple(UE,RIS,signal,y,W,regime)
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    Delta=0.01;     % resolution of 1 cm        
    Pnum=round((UE.rhoRange(2)-UE.rhoRange(1))/Delta);
    Pnum=min(2000,Pnum);   
    rho_grid=linspace(UE.rhoRange(1),UE.rhoRange(2),Pnum);              
    PLogweight=zeros(1,Pnum);       
    % generate possible locations
    Ploc=rho_grid.*[cos(UE.bestPhi).*sin(UE.bestTheta); sin(UE.bestPhi).*sin(UE.bestTheta); cos(UE.bestTheta)];            
    for k=1:Pnum            
        % compute likelihood
        [~,RISphases]=computeRISChannel(Ploc(:,k),RIS,signal,regime);        
        PLogweight(k)=computeLogLikelihood(exp(-1j*RISphases),y,W,signal);           
    end                
    [LLF,index]=max(PLogweight);
    Pest=Ploc(:,index);     
   
