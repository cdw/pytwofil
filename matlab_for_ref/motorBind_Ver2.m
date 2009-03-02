function [Mf,Af] = motorBind_Ver2(Mf, Af, Sc)

% [Mf,Af] = motorBind_Ver1(Mf, Af)
% transitions XB linkages between states
%
% Mf - the myosin filament structure
% Af - the actin filament structure

%% General comments and doc

% CDW(20070829)-This is more or less along the lines of the previous bind state
%               code. This means that it does not actually use proper transition 
%               rates but has this probability method that is not so great. This
%               is ok but has to be fixed before used for a Biophys poster.
% CDW(20080110)-The transition rates are now fixed for the Biophys poster, at
%               least to some extent, they should still be looked at more

%% Code

%Uncomment this next line when a new ver is saved
%warning(['Running an old version of ' mfilename]) 


% set up some constants that we will use in the rest of this function
RVals = rand(1, length(Mf.bst));    %random values to be used in transitions


%for all motors
for m=1:length(Mf.bst),

    %if the motor is unbound
    if Mf.bst(m)==0,
        %find the closest actin site
        [MinVal,MinInd] = min(sqrt(...
        	(Af.loc(1,:)-Mf.hloc(1,m)).^2 + ...
            (Af.loc(2,:)-Mf.hloc(2,m)).^2   ));
        %and if that site is unbound and hurdles over the probability of binding
        if 1-(exp(-MinVal^2)) < RVals(m) && Af.bst(MinInd)==0
            %link the thick fil
            Mf.lnk(:,m) = MinInd;
            Mf.bst(m)   = 1;
            %set thick fil values
            Mf.rs(m)    = Mf.rv(1,2);
            Mf.ths(m)   = Mf.thv(1,2);
            Mf.phs(m)   = Mf.phv(1,2);
            Mf.rk(m)    = Mf.rv(2,2);
            Mf.thk(m)   = Mf.thv(2,2);
            Mf.phk(m)   = Mf.phv(2,2);
            %link the thin fil
            Af.lnk(:,MinInd) = m;
            Af.bst(MinInd)   = 1;
        end

    %or if the motor is loosely bound
    elseif Mf.bst(m)==1,
        %find the distance between the A and M nodes
        Dist = sqrt(...
        	(Af.loc(1,Mf.lnk(m))-Mf.hloc(1,m)).^2 + ...
            (Af.loc(2,Mf.lnk(m))-Mf.hloc(2,m)).^2   );
        %it can bind more tightly
        if 1-.001*100/sqrt(2*0.2515)*(1-tanh(1*sqrt(2*0.2515)*(Dist-6))) < RVals(m) %rand < .1,
            %set thick fil binding state
        	Mf.bst(m) = 2;
            %and thick fil values
            Mf.rs(m)  = Mf.rv(1,3);
            Mf.ths(m) = Mf.thv(1,3);
            Mf.phs(m) = Mf.phv(1,3);
            Mf.rk(m)  = Mf.rv(2,3);
            Mf.thk(m) = Mf.thv(2,3);
            Mf.phk(m) = Mf.phv(2,3);
            %record the force generation as atp useage
            %we may change this in the future to allow recording of all binding
            %state transitions, but for now this will do - CDW (20080114)
            Mf.atp(m) = 1;
        %or unbind
        elseif (exp(15*Mf.rk(m)*0.2515*(...   e^(a*r_k,RT*
                sqrt((Af.loc(1,Mf.lnk(m))-Mf.loc(1,m))^2 + (Af.loc(2,Mf.lnk(m))-Mf.loc(2,m))^2)... (r
                - Mf.rs(m))^2 + ... -rs)^2
                3*(Mf.thk(m)*0.2515)/... + r*th_k,RT /
                sqrt((Af.loc(1,Mf.lnk(m))-Mf.loc(1,m))^2 + (Af.loc(2,Mf.lnk(m))-Mf.loc(2,m))^2)... r
               *atan2((Af.loc(2,Mf.lnk(m))-Mf.loc(2,m)) , (Af.loc(1,Mf.lnk(m))-Mf.loc(1,m)))... *(th-
               - Mf.ths(m))^2 - 1) < RVals(m), %th_s)^2 -1 
            %unlink the thin fil
            Af.bst(Mf.lnk(m)) = 0;
            Af.lnk(Mf.lnk(m)) = 0;
            %reset thick fil values
            Mf.rs(m)  = Mf.rv(1,1);
            Mf.ths(m) = Mf.thv(1,1);
            Mf.phs(m) = Mf.phv(1,1);
            Mf.rk(m)  = Mf.rv(2,1);
            Mf.thk(m) = Mf.thv(2,1);
            Mf.phk(m) = Mf.phv(2,1);
            %unlink the thick fil
            Mf.lnk(m) = 0;
            Mf.bst(m) = 0;
        end
        
    %or if the motor is tightly bound
    elseif Mf.bst(m)==2,
        %find the distance and the angle of sep
        Dist = sqrt((Af.loc(1,Mf.lnk(m))-Mf.loc(1,m))^2 + (Af.loc(2,Mf.lnk(m))-Mf.loc(2,m))^2);
        Ang = atan2((Af.loc(2,Mf.lnk(m))-Mf.loc(2,m))  ,  (Af.loc(1,Mf.lnk(m))-Mf.loc(1,m)));
        %it can bind more loosely
        % if .005 /exp(E2-E3)<random value
        if (0.005* 1/exp(...
                Mf.rv(2,2)*0.2515 *(Dist-Mf.rv(1,2))^2 + ...
                (Mf.thv(2,2)*0.2515/Dist)*(Ang-Mf.thv(1,2))^2 - ...
                Mf.rv(2,3)*0.2515 *(Dist-Mf.rv(1,3))^2 + ... % I think this should be a - at the end CDW20090225
                (Mf.thv(2,3)*0.2515/Dist)*(Ang-Mf.thv(1,3))^2 ...
                )) < RVals(m),
            Mf.bst(m)   = 1;
            %set thick fil values
            Mf.rs(m)    = Mf.rv(1,2);
            Mf.ths(m)   = Mf.thv(1,2);
            Mf.phs(m)   = Mf.phv(1,2);
            Mf.rk(m)    = Mf.rv(2,2);
            Mf.thk(m)   = Mf.thv(2,2);
            Mf.phk(m)   = Mf.phv(2,2);
        %or unbind
        elseif (.001*(sqrt(Mf.rk(m)*0.2515)*(sqrt(3600*Dist^2) - 40*Dist)+ 20)) > RVals(m),
            %unlink the thin fil
            Af.bst(Mf.lnk(m)) = 0;
            Af.lnk(Mf.lnk(m)) = 0;
            %reset thick fil values
            Mf.rs(m)  = Mf.rv(1,1);
            Mf.ths(m) = Mf.thv(1,1);
            Mf.phs(m) = Mf.phv(1,1);
            Mf.rk(m)  = Mf.rv(2,1);
            Mf.thk(m) = Mf.thv(2,1);
            Mf.phk(m) = Mf.phv(2,1);
            %unlink the thick fil
            Mf.lnk(m) = 0;
            Mf.bst(m) = 0;
        end
        
    end %end of if series for each motor


end %end of for loop over all motors
