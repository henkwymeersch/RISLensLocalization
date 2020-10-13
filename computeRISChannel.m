function [gain, phase_rot, locations]=computeRISChannel(source,IRS,signal,regime)
% function [gain, phase_rot, locations]=computeRISChannel(source,RIS,signal,regime)
% computes nearfield of farfield complex channel gain vector
%
% inputs: 
%   source: 3D transmiter location
%   IRS: RIS data structure
%   signal: signal data structure
%   regime: different channel models
%
% outputs:
%   complex gain across the RIS
%   -phases across the RIS
%   RIS element locations
% (c) 2020, Henk Wymeersch, henkw@chalmers.se

M=IRS.M;
fc=signal.fc;
c=signal.c;           
spacing=IRS.Delta;
RIS=IRS.Location;
lambda=c/fc;     % Wavelenght(m)
A=(spacing)^2;   % basic element area
a=sqrt(A);       % basic element size

phi=atan2(source(2),source(1));
theta=acos(source(3)/norm(source));
k=2*pi/lambda*[cos(phi)*sin(theta); sin(phi)*sin(theta) ; cos(theta)]; % wavevector
gain=zeros(M,1);
phase_rot=zeros(M,1);
iix=0:sqrt(M)-1;
iiy=0:sqrt(M)-1;
iix=a*(iix-sqrt(M)/2);
iiy=a*(iiy-sqrt(M)/2);
locations=zeros(3,M);
coord=repmat(1:sqrt(M),sqrt(M),1);
Xtemp=iix(coord);
Ytemp=iiy(coord');
XYZmat=[Xtemp(:)';Ytemp(:)'];
locations(1:2,:)=XYZmat;
d=vecnorm(locations-source);
d0=norm(source);
% Correction factor in order to match the nearfield power. 
correction=1-sin(theta)^2*sin(phi)^2;
switch regime
    case 'CM1'         
        phase_rot=mod(-k'*locations,2*pi);
        phase_rot=phase_rot';
        gain=((sqrt(cos(theta)*correction)*a)/(sqrt(4*pi)*norm(source-RIS)))*exp(-1i*phase_rot);
    case 'CM2'
        phase_rot=mod(2*pi*(d-d0)/lambda,2*pi);        
        phase_rot=phase_rot';
        gain=((sqrt(cos(theta)*correction)*a)/(sqrt(4*pi)*norm(source-RIS)))*exp(-1i*phase_rot);
    case 'CM3'        
    for m=1:M          
        el_m =locations(:,m);
        x_m=el_m(1);
        y_m=el_m(2);       
        setX=[a/2+x_m-source(1) a/2-x_m+source(1)];             
        setY=[a/2+y_m-source(2) a/2-y_m+source(2)];                                              
        d=abs(source(3));                            % according to the paper, the Z-coordinate is d
        TMP=zeros(length(setX),length(setY));
        dm=norm(source-el_m);                   % distance between BS and m-th element
        for ix=1:length(setX)
            for iy=1:length(setY)                        
                a2=(setX(ix)*setY(iy))/d^2;                        
                b21=3*((setY(iy)^2)/d^2)+3;                        
                b22=sqrt(((setX(ix)^2+setY(iy)^2)/d^2)+1);                       
                TMP(ix,iy)=(a2/(b21*b22))+((2/3)*atan2(a2,b22));                        
            end
        end                     
        powerval=((1/(4*pi))*sum(sum(TMP)));                % power of the m-th element                         
        phase_rot(m)=2*pi*(dm-d0)/lambda;
        gain(m)=sqrt(powerval)*exp(-1i*phase_rot(m));     % complex channel gain                               
    end
end