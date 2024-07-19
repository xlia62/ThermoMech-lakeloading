function MECL = eclogite(RxType,MTK,MPR,MI,MECL,MPHASE,RGAS,timestep)

% Function to model instantaneous phase transition from basalt to eclogite
% We follow the approach outlined in Giunchi and Ricard (GJI 1998) and 
% van Hunen et al. (Tectonophysics, 2002) approach without solving the kinetics
%
% albite <-> jadeite + quartz
%
% MECL = eclogite(RxType,MTK,MPR,MI,MECL,MPHASE,timestep)
%
% MTK           = marker temperature
% MPR           = marker total pressure
% MI            = marker id rock type
% MECL(:,1)     = Fraction of eclogite
% MECL(:,2)     = Previous change in Gibb's free energy
% MECL(:,3)     = Previous reaction rate
% MECL(:,4)     = Reaction time
% MECL(:,5)     = Starting fraction of basalt (prograde) or eclogite (retrograde)
% MPHASE(1,:) = Kinetic Activation energy [prograde retrograde]
% MPHASE(2,:) = Kinetic Pre-exponential factor [prograde retrograde]
% MPHASE(3)   = Avrami exponent,
% MPHASE(4)   = Non kinetic phase transformation Pressure
%               = 4 means rate of transformation is not constant, but starts 
%                 at zero and ends at zero
%               = 1 reaction is fast at first than slow
% timestep      = time step size
%
% V1.0 Robert Moucha, Syracuse University, February 25, 2017
%

switch RxType

    case {'Depth','depth'}  % Prograde only
        
        % Transition occures at provided transition pressure Ptr (i.e. depth)       
        
        k = (MI == 4 | MI == 5 | MI == 7) & MPR > MPHASE(4,1);
        
        MECL(k,1) = 1;          % Fraction of transformation = 100%
                
    case {'Thermodynamic','thermodynamic'}  % Prograde only
        
        % Transition occures when the equilibrium of the eclogite phase (albite=jadeite+quartz reaction) is
        % crossed combined with a phase boundary from Spears 1993
        % depends on both Temperature and Pressure       
        
        k = find((MI == 4 | MI == 5 | MI == 7) & MECL(:,1) < 1);
        
        dG = eclogiteGfree(MTK(k),MPR(k));
        
        % Stable eclogite field
        j = dG < 0;
        
        % Now determine if the phase boundary from blueschist to eclogite
        % was crossed.
        
        % Blueschist to eclogite boundary (linear approx., GPa and oC) according to Spear
        % 1993
        p(1) = -0.0030 ;
        p(2) = 2.26;
        
        pr = p(1)*(MTK(k)-273)+p(2);
        j2 = MPR(k)/1e9 > pr & j;
        
        MECL(k(j2),1) = 1;          % Fraction of transformation = 100%
        
        % Future things to include
        
    case {'Kinetic','kinetic'}
        
        % Use kinetic reaction rate to compute fraction of eclogite
        k = find(MI == 4 | MI == 5 | MI == 7);
        nk = length(k);
        
        T = MTK(k);         % Temperature
        P = MPR(k);         % Pressure
        F1 = MECL(k,1);     % Fraction of eclogite
        dG1 = MECL(k,2);    % Previous change in Gibb's free energy
        Y1 = MECL(k,3);     % Previous reaction rate
        t1 = MECL(k,4);     % Reaction time
        F0 = MECL(k,5);     % Starting fraction of basalt (prograde) or eclogite (retrograde)
        
        % Set to retrograde default kinetics, use a vector.
        Ek = MPHASE(1,2) + zeros(nk,1);  
        Ak = MPHASE(2,2) + zeros(nk,1);   
        n = MPHASE(3,1);  % Same for both directions
        
        % Compute the change in Gibb's free energy
        dG = eclogiteGfree(T,P);
          
        % Is there a change in direction of phase transformation? 
        change_dir = sign(dG) ~= sign(dG1);
        
        % If change in Gibbs free energy is negative it is prograde = true
        % (retrograde is false)
        prograde = dG < 0;
        Ek(prograde) = MPHASE(1,1);
        Ak(prograde) = MPHASE(2,1);

        % Reaction rate
        Y2 = T.*sign(-dG).*exp(-Ek./(RGAS*T)).*(1-exp(abs(dG)./(RGAS*T)));
        
        % Compute reaction time
        t = (Y1./Y2).*t1 + timestep;
        % Reset reaction time when direction changes
        t(change_dir) = timestep;
        
        % Reset the concentrations if there is a change in direction
        F0(change_dir) = F1(change_dir);
        F0(change_dir & prograde) = 1 - F1(change_dir & prograde);
        
        % Compute the cummulative fraction for retrograde
        F = exp(-(Ak.*Y2.*t).^n);
        
        % If prograde change the fraction
        F(prograde) = 1 - F(prograde);
        
        % Account for initial concentration (fraction)
        F = F0.*F;
        
        % Put back into container for output
        MECL(k,1) = F;
        MECL(k,2) = dG;
        MECL(k,3) = Y1;
        MECL(k,4) = t;
        MECL(k,5) = F0;
        
    otherwise
        error('Unknown Phase Transition Reaction Type')
end
