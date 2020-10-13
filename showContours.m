function showContours(signal,UE,RIS)
% PEB computation in space
% (c) 2020, Henk Wymeersch, henkw@chalmers.se

    x_grid=linspace(-2,2,50);
    z_grid=linspace(0.001,4,50);
    Omega_b=diag(conj(RIS.b./(abs(RIS.b))));
    bnew=RIS.b.'*Omega_b;    
    PEB=zeros(length(x_grid),length(z_grid));
    SNRr=zeros(length(x_grid),length(z_grid));
    SNRd=zeros(length(x_grid),length(z_grid));
    SNRp=zeros(length(x_grid),length(z_grid));

    UE1=UE;
    UE1.Location=0.1*[1;1;1];       % assume beams are pointing here        
    UE1.rho=norm(UE1.Location);
    UE1.phi=atan2(UE1.Location(2),UE1.Location(1));  % between 0 and 2pi
    UE1.theta=acos(UE1.Location(3)/norm(UE1.Location)); % between 0 and pi (for us pi/2, since Z>0)  
    UE1.mean=UE1.Location;   
    UE1.covariance=UE1.covariance;        
    BeamsR=getBeams(UE1,RIS,signal,'random');
    BeamsD=getBeams(UE1,RIS,signal,'direction');        
    BeamsP=getBeams(UE1,RIS,signal,'position');        
    UE.Location=zeros(3,1);
  
    for k=1:length(x_grid)                  
        UE.Location(1)=x_grid(k);
        UE.Location(2)=x_grid(k);
        for l=1:length(z_grid)        
            UE.Location(3)=z_grid(l);
            UE.rho=norm(UE.Location);
            UE.phi=atan2(UE.Location(2),UE.Location(1));  % between 0 and 2pi
            UE.theta=acos(UE.Location(3)/norm(UE.Location)); % between 0 and pi (for us pi/2, since Z>0)       
            h=computeRISChannel(UE.Location,RIS,signal,'CM3');           
                        
            % random beams
            W=bnew.*BeamsR;                                 
            SNRr(k,l)=10*log10(norm(sqrt(signal.P)*W*h)^2/(2*signal.sigma2*signal.T));                                   
            
            % directional beams
            W=bnew.*BeamsD;                           
            SNRd(k,l)=10*log10(norm(sqrt(signal.P)*W*h)^2/(2*signal.sigma2*signal.T));            
            
            % positional beams
            W=bnew.*BeamsP;                          
            SNRp(k,l)=10*log10(norm(sqrt(signal.P)*W*h)^2/(2*signal.sigma2*signal.T));            

        end
    end
  
        figure(2);
        subplot(1,3,1)
        SNRth=0;    % SNR threshold
        f1=contourf(x_grid,z_grid,max(SNRr',SNRth),'edgecolor','none');
        xl=xlabel('$x$ [m]');
        yl=ylabel('$z$ [m]');
        tt=title('random, SNR [dB]');
        set(xl,'Interpreter','latex','FontSize',12);
        set(yl,'Interpreter','latex','FontSize',12);
        set(tt,'Interpreter','latex','FontSize',12);
        pbaspect([2 2 1])
        c = colorbar;
        

        subplot(1,3,2)
        f2=contourf(x_grid,z_grid,max(SNRd',SNRth),'edgecolor','none');
        xl=xlabel('$x$ [m]');
        yl=ylabel('$z$ [m]');
        tt=title('directional, SNR [dB]');
        set(xl,'Interpreter','latex','FontSize',12);
        set(yl,'Interpreter','latex','FontSize',12);
        set(tt,'Interpreter','latex','FontSize',12);
        pbaspect([2 2 1])
        c = colorbar;
        


        subplot(1,3,3)
        f3=contourf(x_grid,z_grid,max(SNRp',SNRth),'edgecolor','none');
        xl=xlabel('$x$ [m]');
        yl=ylabel('$z$ [m]');
        tt=title('positional, SNR [dB]');
        set(xl,'Interpreter','latex','FontSize',12);
        set(yl,'Interpreter','latex','FontSize',12);
        set(tt,'Interpreter','latex','FontSize',12);
        c = colorbar;
        pbaspect([2 2 1])
        set(gcf, 'Color', 'w');
    
    
    
