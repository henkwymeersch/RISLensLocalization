function showPEB(signal,UE,RIS,e_unit)
% PEB in line
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
signal.sims=2;
PEBr=zeros(length(UE.distances),signal.sims);
PEBd=zeros(length(UE.distances),signal.sims);
PEBp=zeros(length(UE.distances),signal.sims);
UE.covariance=(eye(3));
for UEDi=1:length(UE.distances)        
     for sim=1:signal.sims
        dd=UE.distances(UEDi);
        % UE parameters 
        % --------------
        UE.Location=dd*e_unit; 
        UE.rho=norm(UE.Location);
        UE.phi=atan2(UE.Location(2),UE.Location(1));  % between 0 and 2pi
        UE.theta=acos(UE.Location(3)/norm(UE.Location)); % between 0 and pi (for us pi/2, since Z>0)
        UE.mean=UE.Location;          
        % Generate observation
        % --------------------
        h=computeRISChannel(UE.Location,RIS,signal,'CM2');
        Omega_b=diag(conj(RIS.b./(abs(RIS.b))));
        bnew=RIS.b.'*Omega_b;
        
        Beams=getBeams(UE,RIS,signal,'random');
        W=bnew.*Beams;
        PEBr(UEDi,sim)=computePEB(UE,RIS,signal,W,'CM2');         
        
         Beams=getBeams(UE,RIS,signal,'direction');
        W=bnew.*Beams;            
        PEBd(UEDi,sim)=computePEB(UE,RIS,signal,W,'CM2');         
        
        Beams=getBeams(UE,RIS,signal,'position');
        W=bnew.*Beams;                   
        PEBp(UEDi,sim)=computePEB(UE,RIS,signal,W,'CM2'); 
        
        
        % now also with concentrated a priori information
        UE1=UE;
        UE1.covariance=0.01*(eye(3));
        
        Beams=getBeams(UE1,RIS,signal,'direction');
        W=bnew.*Beams;            
        PEBd2(UEDi,sim)=computePEB(UE1,RIS,signal,W,'CM2');         
        
        Beams=getBeams(UE1,RIS,signal,'position');
        W=bnew.*Beams;                   
        PEBp2(UEDi,sim)=computePEB(UE1,RIS,signal,W,'CM2'); 

        
     end
end

figure(1)
semilogy(UE.distances,mean(PEBr'),'s-',UE.distances,mean(PEBd'),'+-',UE.distances,mean(PEBp'),'*-',UE.distances,mean(PEBd2'),'+--',UE.distances,mean(PEBp2'),'*--')
xlabel('distance [m]')
ylabel('PEB [m]')
grid
legend('randomized','directional, \sigma = 1 ','positional, \sigma = 1','directional, \sigma = 0.1','positional, \sigma = 0.1','FontSize',12)


