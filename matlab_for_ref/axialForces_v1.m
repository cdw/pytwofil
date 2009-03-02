function[MfForces]=axialForces_v1(Mf, Af, Sc)
%the purpose of this function is to calculate the axial forces experienced by
%each of the thick filament nodes, this is all ripped from balFils_v1

c = [Mf.loc(1,:),Af.loc(1,:)];   %setup c
Mn = length(Mf.loc);
An = length(Af.loc);
Sep = Sc.sep;
Len = Sc.len;


    %force for first point of the thick filament
    MfForces(1)=Mf.mk*(c(2)-c(1)-Mf.ms)- Mf.mk*(c(1)-0-Mf.ms)  ... % force from adjacent site
        + sign(Mf.bst(1))*... %turn on/off the XB force as a function of the XB state
        (Mf.rk(1)* (cLenXbMf(1)-Mf.rs(1))* cos(cThMf(1)) ... %Kxb force
        - (1/cLenXbMf(1))* Mf.thk(1)* (cThMf(1)-Mf.ths(1))* sin(cThMf(1))); %Kth force

    %force on each point ('cept first and last)
    MfForces(2:Mn-1)=Mf.mk.*(c(3:Mn)-c(2:Mn-1)-Mf.ms) - Mf.mk.*(c(2:Mn-1)-c(1:Mn-2)-Mf.ms) ... %force of adjacent sites
        	+ sign(Mf.bst(2:Mn-1)).*... %turn on/off the XB force as a function of the XB state
            (Mf.rk(2:Mn-1).* (cLenXbMf(2:Mn-1)-Mf.rs(2:Mn-1)).* cos(cThMf(2:Mn-1)) ... %Kxb force
            - (1./cLenXbMf(2:Mn-1)).* Mf.thk(2:Mn-1).* (cThMf(2:Mn-1)-Mf.ths(2:Mn-1)).* sin(cThMf(2:Mn-1))); %Kth force

    %calculate the force at the last point of the thick filament
    MfForces(Mn)=-Mf.mk*(c(Mn)-c(Mn-1)-Mf.ms) ... % force from adjacent thick sites
        + sign(Mf.bst(Mn))*... %moderate xb force with xb state
        (Mf.rk(Mn)* (cLenXbMf(Mn)-Mf.rs(Mn))* cos(cThMf(Mn)) ... %Kxb force
        - (1/cLenXbMf(Mn))* Mf.thk(Mn)* (cThMf(Mn)-Mf.ths(Mn))* sin(cThMf(Mn))); %Kth force
    




%% Subsidiary functions to make stuff easier to implement above

    function [Theta] = cThMf(index)
    %this calculates the angle theta at any Mf node, and
    %returns zero if that index is not linked to any other
        Theta = atan2(Sep, ...    %take the arctangent of the separation and x
            (c(Mn+Mf.lnk(index))-c(index)));  
    end %ends cTh

%     function [Theta] = cThAf(index)
%     %this calculates the angle theta at any Af node, and
%     %returns zero if that index is not linked to any other
%         Theta = atan2(Sep, ...    %take the arctangent of the separation and x
%             (c(index)-c(chZ(Af.lnk(index-Mn)))));  %or the thin fil
%     end %ends cTh


    function [Length] = cLenXbMf(index)
    Length = ...(index<=Mn)... if we are looking at the thick fils      	(sign(c(Mn+Mf.lnk(index))-c(index))+((c(Mn+Mf.lnk(index))-c(index))==0)).* ... have we crossed over a node?
         realsqrt(realpow(c(Mn+Mf.lnk(index))-c(index),2)+realpow(Sep,2));%... what is the distance
    end %ends cLenXB

%     %calculates the length of any xb on af, returning -1 if there is no link
%     function [Length] = cLenXbAf(index)
%         Length = ...        	(sign(c(index)-c(chZ(Af.lnk(index-Mn))))+((c(index)-c(chZ(Af.lnk(index-Mn))))==0)).* ... have we crossed over a node?
%             realsqrt(realpow(c(index)-c(chZ(Af.lnk(index-Mn))),2)+realpow(Sep,2)); % and what is the distance
%     end %ends cLenXB




end