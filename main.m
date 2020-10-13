% --------------------------------------------
% MATLAB code related to the paper
% Near-field Localization with a Reconfigurable Intelligent Surface Acting as Lens
% https://arxiv.org/abs/2010.05617
% Zohair Abu-Shaban, Kamran Keykhosravi, Musa Furkan Keskin, George C. Alexandropoulos, Gonzalo Seco-Granados, Henk Wymeersch
% code, (c) 2020, Henk Wymeersch, henkw@chalmers.se
% --------------------------------------------


close all
clear all
warning off
clc


% 1. signal parameters
% ---------------------
signal.fc=28;                       % carrier in GHz
signal.c=0.3;                       % speed of light [m/ns]
signal.T=200;                       % number of observations
signal.P=1;                         % transmit power mW
signal.N0dBmHz=-174;                % dBm/Hz
signal.N0=10^(0.1*signal.N0dBmHz)*1e9; % noise PSD  [mW/GHz] (290 Kelvin * Boltzmann constant in W/Hz)
signal.BW=1e-3;                     % Bandwidth GHz
signal.NFdB=8;                      % receiver noise figure [dB]
signal.NF=10^(0.1*signal.NFdB);
signal.EsN0=2*signal.P/(signal.NF*signal.N0*signal.BW);
signal.sigma2=signal.NF*signal.N0*signal.BW/2;
signal.lambda=signal.c/signal.fc;   % wavelength
signal.plot=0;                      % 1 = show plots
signal.sims=5;                      % number of Monte Carlo simulations

% 2. RIS parameters
% --------------
RIS.M=50*50;                            % RIS elements 
RIS.Delta=signal.lambda/2;              % RIS element spacing
RIS.Location=[0;0;0];                   % location of RIS in XY plane
RIS.Antenna=[0; 0; -1*signal.lambda];   % RIS antenna location
RIS.Orientation=0;                      % RIS rotation around the Y axis 
RIS.b=computeRISChannel(RIS.Antenna,RIS,signal,'CM3');   % response from RIS to antenna
RIS.bits=1000;
RIS.beamType='random'; % 'position', 'direction', 'random'



% 3. UE parameters 
% --------------
UE.covariance=0.01*diag([1;1;1]);   % UE a priori covariance
UE.distances=[0.15 1 1.5 2 3 4 6 8 10 12 15];
e_unit=[1; 1; 1]/norm([1; 1; 1]);


% these two function show the first and second figure from the paper
showContours(signal,UE,RIS)
showPEB(signal,UE,RIS,e_unit)


% wipe of phase of b
ab=quantizePhases(angle(RIS.b),RIS.bits);        
Omega_b=diag(exp(-1j*ab));        
bnew=RIS.b.'*Omega_b;
errorStat=zeros(length(UE.distances),signal.sims);
RMSE=zeros(1,length(UE.distances));
for UEDi=1:length(UE.distances)        
    dd=UE.distances(UEDi);
    for sim=1:signal.sims        
        text=['simulation ' num2str(sim) ' for distance ' num2str(dd)];
        disp(text)
        % UE parameters 
        % --------------
        UE.Location=dd*e_unit; 
        UE.rho=norm(UE.Location);
        UE.phi=atan2(UE.Location(2),UE.Location(1));        % between 0 and 2pi
        UE.theta=acos(UE.Location(3)/norm(UE.Location));    % between 0 and pi (for us pi/2, since Z>0)        
        UE.mean=mvnrnd(UE.Location,UE.covariance)';           
        UE.mean(3)=abs(UE.mean(3)); % Z coordinate is positive
        
        % Generate observation
        % --------------------
        h=computeRISChannel(UE.Location,RIS,signal,'CM3');                                
        Beams=getBeams(UE,RIS,signal,RIS.beamType);
        W=bnew.*Beams;
        noise=randn(signal.T,1)+1j*randn(signal.T,1);
        s=sqrt(signal.P)*W*h;
        y=s+noise*sqrt(signal.sigma2);  
        
        % Estimate the position
        % ----------------------
        UE.rhoRange=[max(0,UE.rho-5) UE.rho+5];                
        UE.phiRange=[0 2*pi];
        UE.thetaRange=[0 pi/2];        
        UE=getAngleEstimateSimple(UE,RIS,signal,y,W);                               
        [UE.estimate,LLF2]=getPositionSimple(UE,RIS,signal,y,W,'CM2');        
        errorStat(UEDi,sim)=norm(UE.estimate-UE.Location);               
    end
    yr=errorStat(UEDi,:).^2;    
    RMSE(UEDi)=sqrt(mean(yr'));  
    figure(5)
    semilogy(UE.distances,RMSE)
    ylabel('RMSE')
    xlabel('distance to RIS')
    grid
end

