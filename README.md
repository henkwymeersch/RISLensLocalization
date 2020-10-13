# RIS Lens Localization
3D Localization with RIS lens and a single antenna

## Summary
The matlab code (main.m) generates a 3D environment with an intelligent reconfigurable lens and a transmitter. The transmitter sends a narrowband signal to the lens. Different phase profiles at the lens are used to generate a scalar stream of observations, from which the transmitter location is estimated. The script generates 3 figures from the following paper: 
Zohair Abu-Shaban, Kamran Keykhosravi, Musa Furkan Keskin, George C. Alexandropoulos, Gonzalo Seco-Granados, Henk Wymeersch, "Near-field Localization with a Reconfigurable Intelligent Surface Acting as Lens", arXiv:2010.05617.

## Main parameters
The parameters are the following
```
% 1. signal parameters
% ---------------------
signal.fc=28;                       % carrier in GHz
signal.c=0.3;                       % speed of light [m/ns]
signal.T=200;                       % number of observations
signal.P=1;                         % transmit power mW
signal.N0dBmHz=-174;                % dBm/Hz
signal.BW=1e-3;                     % Bandwidth GHz
signal.NFdB=8;                      % receiver noise figure [dB]

% 2. RIS parameters
% --------------
RIS.M=50*50;                            % RIS elements 
RIS.Delta=signal.lambda/2;              % RIS element spacing
RIS.Location=[0;0;0];                   % location of RIS in XY plane
RIS.Antenna=[0; 0; -1*signal.lambda];   % RIS antenna location
RIS.bits=1000;                          % number of bits per RIS element
RIS.beamType='random'; % 'position', 'direction', 'random'

% 3. UE parameters 
% --------------
UE.covariance=0.01*diag([1;1;1]);   % UE a priori covariance
UE.distances=[0.15 1 1.5 2 3 4 6 8 10 12 15]; % UE distances
```
These paramaters can be changed to obtain different results. 
