function [Mf, Af, exitflag] = balFils_v4(Mf,Af,Sc,Method)

%[Mf, Af, exitflag] = balFils_v2(Mf,Af,Sc,Method)
%
%This function serves to perform a force balance on the two fil system
% Ver_4 is like Ver_2, but with a rewritten force equation

%% Old verision protection
%Uncomment this next line when a major new ver is saved
% if ~strcmp(lastwarn, ['Running an old version of ' mfilename]) %only warn once
%   warning('DangerDave:OldVersion',['Running an old version of ' mfilename]) 
% end


%% Setup and inititate the fsolve equation

c = [Mf.loc(1,:),Af.loc(1,:)];   %setup c
Mn = length(Mf.loc);
An = length(Af.loc);
Sep = Sc.sep;
Len = Sc.len;

%unpack more variables and matrices in an attempt to imporve readability
Mk  = Mf.mk;
Ms  = Mf.ms;
Uda = Mf.uda;
Mln = Mf.lnk;
Mbs = Mf.bst;
Rk  = Mf.rk;
Rs  = Mf.rs;
THk = Mf.thk;
THs = Mf.ths;
Ak  = Af.ak;
As  = Af.as;
Aln = Af.lnk;
Abs = Af.bst; 

%we are commenteing out the system we actually want to use- 6-12-07
options = optimset('NonlEqnAlgorithm', Method,...   %using the gauss newton alg                'Diagnostics','on',...
                   'Display','none',...
                   'FunValCheck','on',...                   'PlotFcns','@optimplotfval',...
                   'TolFun',1e-5);%,...                   'PlotFcns' ,{  @optimplotx @optimplotfval });
               
[c,fval,exitflag] = fsolve(@forceForFSolve,c,options); %solve for root


%update the filiament locations
Mf.loc(1,:) = c(1:Mn);
Af.loc(1,:) = c(Mn+1:Mn+An);

%% Start of force calculating equations

function F = forceForFSolve(c)
    % Shoot, gonna have to give a good explanation of c here,
    % c is all the x locations of the thick and thin filaments e.g.
    % c = [ThickXLoc1, ThickXLoc2,...,ThinXLoc1,ThinXLoc2...]

    %Preallocate the matrix of output forces
    F = zeros(Mn+An,1);
    
%% Forces on thick filament sites

    %force for first point of the thick filament
    if Mf.bst(1)==0
        F(1)=Mk*(c(2)-c(1)-Ms)- Mk*(c(1)-0-Uda); % force from adjacent site
    else
        MLenXb = realsqrt(realpow(c(Mn+Mln(1))-c(1),2)+realpow(Sep,2));
        MTh    = atan2(Sep, (c(Mn+Mln(1))-c(1)));
        F(1)=Mk*(c(2)-c(1)-Ms)- Mk*(c(1)-0-Uda) ... % force from adjacent site
            + Rk(1)* (MLenXb-Rs(1))* cos(MTh) ... %Kxb force
            - (1/MLenXb)* THk(1)* (MTh-THs(1))* sin(MTh); %Kth force
    end

    %force on each point ('cept first and last)
    for i = 2:Mn-1
        if Mf.bst(i)==0
            F(i)=Mk*(c(i+1)-c(i)-Ms) - Mk*(c(i)-c(i-1)-Ms); %force of adjacent sites
        else
            MLenXb = realsqrt(realpow(c(Mn+Mln(i))-c(i),2)+realpow(Sep,2));
            MTh    = atan2(Sep, (c(Mn+Mln(i))-c(i)));
            F(i)=Mk.*(c(i+1)-c(i)-Ms) - Mk.*(c(i)-c(i-1)-Ms) ... %force of adjacent sites
                + Rk(i).* (MLenXb-Rs(i)).* cos(MTh) ... %Kxb force
                - (1./MLenXb).* THk(i).* (MTh-THs(i)).* sin(MTh); %Kth force
        end
    end

    %calculate the force at the last point of the thick filament
    if Mf.bst(20)==0
        F(Mn)=-Mk*(c(Mn)-c(Mn-1)-Ms); % force from adjacent thick sites
    else
        MLenXb = realsqrt(realpow(c(Mn+Mln(Mn))-c(Mn),2)+realpow(Sep,2));
        MTh    = atan2(Sep, (c(Mn+Mln(Mn))-c(Mn)));
        F(Mn)=-Mk*(c(Mn)-c(Mn-1)-Ms) ... % force from adjacent thick sites
            + Rk(Mn)* (MLenXb-Rs(Mn))* cos(MTh) ... %Kxb force
            - (1/MLenXb)* THk(Mn)* (MTh-THs(Mn))* sin(MTh); %Kth force
    end
    
%% Forces on thin filament sites
    %calculate the force on the thin filament's first point
    if Af.bst(1)==0
        F(Mn+1)=Ak*(c(Mn+1+1)-c(Mn+1)-As); %<- force from adjacent sites
    else
        ALenXb = realsqrt(realpow(c(Mn+1)-c(Aln(1)),2)+realpow(Sep,2));
        ATh  = atan2(Sep, (c(Mn+1)-c(Aln(1))));
        F(Mn+1)=Ak*(c(Mn+1+1)-c(Mn+1)-As)  ... %<- force from adjacent sites
            - Rk(Aln(1))* (ALenXb-Rs(Aln(1)))* cos(ATh) ... %Kxb force
            + (1/ALenXb)* THk(Aln(1))* (ATh-THs(Aln(1)))* sin(ATh); %Kth force
    end
    
    %calculate force on all points ('cept first and last) of thin filament
    for i = Mn+2:An+Mn-1
        if Af.bst(i-Mn)==0
            F(i)=Ak*(c(i+1)-c(i)-As) - Ak*(c(i)-c(i-1)-As); %force of adj thin sites
        else
            ALenXb = realsqrt(realpow(c(i)-c(Aln(i-Mn)),2)+realpow(Sep,2));
            ATh  = atan2(Sep, (c(i)-c(Aln(i-Mn))));
            F(i)=Ak.*(c(i+1)-c(i)-As) - Ak.*(c(i)-c(i-1)-As) ... %force of adj thin sites
                - Rk(Aln(i-Mn)).*(ALenXb-Rs(Aln(i-Mn))).*cos(ATh) ... %Kxb force
                + (1./ALenXb).* THk(Aln(i-Mn)).* (ATh-THs(Aln(i-Mn))).* sin(ATh); %Kth force
        end
    end

    %calculate the force at the last point of the thin filament
    if Af.bst(30)==0
        F(An+Mn)=Ak*(Len-c(An+Mn)-As)-Ak*(c(An+Mn)-c(An+Mn-1)-As); %force of adjacent thin sites
    else
      ALenXb = realsqrt(realpow(c(An+Mn)-c(Aln(An)),2)+realpow(Sep,2));
      ATh  = atan2(Sep, (c(An+Mn)-c(Aln(An))));
      F(An+Mn)=Ak*(Len-c(An+Mn)-As)-Ak*(c(An+Mn)-c(An+Mn-1)-As) ... %force of adjacent thin sites
            - Rk(Aln(An))* (ALenXb-Rs(Aln(An)))* cos(ATh) ... %Kxb force
            + (1/ALenXb)* THk(Aln(An))* (ATh-THs(Aln(An)))* sin(ATh); %Kth force
    end

    
% %% calculate the xb angle in forceForFSolve
%     function [Theta] = cThMf(index)
%     %this calculates the angle theta at any Mf node, and
%     %returns zero if that index is not linked to any other
%         Theta = atan2(Sep, ...    %take the arctangent of the separation and x
%             (c(Mn+Mln(index))-c(index)));  
%     end %ends cTh
%     
% %% calculate the xb angle in forceForFSolve
%     function [Theta] = cThAf(index)
%     %this calculates the angle theta at any Af node, and
%     %returns zero if that index is not linked to any other
%         Theta = atan2(Sep, ...    %take the arctangent of the separation and x
%             (c(index)-c(Aln(index-Mn))));  %or the thin fil
%     end %ends cTh
% 
%  
% %% calculate the xb length in forceForFSolve
%     %calculates the length of any xb, returning -1 if there is no link
%     function [Length] = cLenXbMf(index)
%     Length = realsqrt(realpow(c(Mn+Mln(index))-c(index),2)+realpow(Sep,2));  %what is the distance
%     end %ends cLenXB
%     
% %% calculate the xb length in forceForFSolve
%     %calculates the length of any xb on af, returning -1 if there is no link
%     function [Length] = cLenXbAf(index)
%         Length = realsqrt(realpow(c(index)-c(Aln(index-Mn)),2)+realpow(Sep,2)); %what is the distance
%     end %ends cLenXB

    
end%ends the calcForceForFSolve

end%ends the balance function
