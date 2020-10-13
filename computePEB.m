function PEB = computePEB(UE,RIS,signal,W,regime)
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    [h,RISphases,locations]=computeRISChannel(UE.Location,RIS,signal,regime);
    a=exp(-1j*RISphases);       
    K=locations-UE.Location;
    d=vecnorm(K);
    K=K./d;
    er=UE.Location/norm(UE.Location);
    Da=1j*2*pi/signal.lambda*(diag(a)*K'+a*er');     
    rho=abs(h(1));
    J=zeros(5,5);   
    for t=1:signal.T                
        myGradient=sqrt(signal.P)*W(t,:)*[a 1j*rho*a rho*Da];
        J=J+real(myGradient'*myGradient);
    end        
    J=J/signal.sigma2;
    Jinv=inv(J);
    PEB=sqrt(trace(Jinv(3:5,3:5)));  