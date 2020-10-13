function UE=getAngleEstimateSimple(UE,RIS,signal,y,W)        
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    % create a grid
    phi_grid=linspace(UE.phiRange(1),UE.phiRange(2),360);
    theta_grid=linspace(UE.thetaRange(1),UE.thetaRange(2),90);                        
    
    % first estimate theta
    [~,~,RISlocations]=computeRISChannel([1;1;1],RIS,signal,'flat');    
    r=vecnorm(RISlocations);
    psi=atan2(RISlocations(2,:),RISlocations(1,:));       
    N=5;        % number of terms in the expansion
    A=zeros(signal.T,2*N+1,length(theta_grid));    
    for k2=1:length(theta_grid)
        theta=theta_grid(k2);
        % make G(theta) matrix
        G=zeros(2*N+1,RIS.M);
        for n=-N:N
            % make 1 row of G
            g=1j^n*exp(-1j*n*psi).*besselj(n,-2*pi*r*sin(theta)/signal.lambda);
            G(n+N+1,:)=g;
        end 
        if (signal.T>N)        
            A(:,:,k2)=W*G.';  % TxN
            B=A(:,:,k2)'*A(:,:,k2); % NxN
            cost(k2)=norm(y-A(:,:,k2)*inv(B)*A(:,:,k2)'*y);
        else
            disp('error: too few observations')
            keyboard
        end
    end
    [mv,mi]=min(cost);
    UE.bestTheta=theta_grid(mi);
    
    % now estimate phi:      
    LogWeight=zeros(1,length(phi_grid));
    for k=1:length(phi_grid)
        p=[sin(UE.bestTheta)*cos(phi_grid(k)); sin(UE.bestTheta)*sin(phi_grid(k)); cos(UE.bestTheta)];                  
        [~,RISphases]=computeRISChannel(p,RIS,signal,'flat');        
        LogWeight(k)=computeLogLikelihood(exp(-1j*RISphases),y,W,signal);       
    end    
    [mv,mi]=max(LogWeight);
    UE.bestPhi=phi_grid(mi);
        