function Beams=getBeams(UE,RIS,signal,beamType)
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    switch beamType
        case 'random'
            phases=rand(signal.T,RIS.M)*2*pi;
            phases=quantizePhases(phases,RIS.bits);            
            Beams=exp(1j*phases);            
        case 'direction'     
            P=mvnrnd(UE.mean,UE.covariance,signal.T)';
            ii=P(3,:)<0;
            P(3,ii)=P(3,ii)*(-1);     % only positive Z axis                           
            for t=1:signal.T            
                [~,RISphases]=computeRISChannel(P(:,t),RIS,signal,'CM1');
                RISphases=quantizePhases(RISphases,RIS.bits);
                a=exp(-1j*RISphases);        
                a=a';       % make conjugate
                Beams(t,:)=a;
            end              
        case 'position'
            P=mvnrnd(UE.mean,UE.covariance,signal.T)';
            ii=P(3,:)<0;
            P(3,ii)=P(3,ii)*(-1);     % only positive Z axis
            for t=1:signal.T
                [~,RISphases]=computeRISChannel(P(:,t),RIS,signal,'CM2');
                RISphases=quantizePhases(RISphases,RIS.bits);
                a=exp(-1j*RISphases);        
                a=a';   % make conjugate
                Beams(t,:)=a;
            end                         
    end