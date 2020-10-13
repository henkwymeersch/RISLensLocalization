function phasesOut=quantizePhases(phasesIn,Nbits)   
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
    if (Nbits>5) % assume perfect resolution
        phasesOut=phasesIn;
    else
        Delta=pi/Nbits;
        phasesOut=floor(phasesIn/Delta+0.5)*Delta;
    end
    