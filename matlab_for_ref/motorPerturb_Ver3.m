function[Mf]=motorPerturb_Ver3(Mf, Sc)

% [Mf]=motorPerturb_Ver1(Mf, Sc)
%
% This finds all unbound heads on a given filament, perturbs them with
% normal random x and y forces, and returns the updated filament with
% the heads in their new positions
% 
% Thick     - A structure containing all info on the thick filament
% ForceMean - The mean of the normal force distribution
% ForceVar  - The variance of the force distribution
% Kxb,Kth   - Spring constants of crossbridge and torsional spring
% Rs, Ths   - Rest length of xbridge and rest angle of spring


%% General Documentation
% CDW(20070731)-This is one of the retrofits from the rotation code
% CDW(20070813)-This is not yet verified to produce proper diffusional
%               forces/energies
% CDW(20080109)-This is now being converted to use more realistic versions of
%               the energies of diffusion in preperation for BPS2008
% CDW(20080110)-The positions throw off imaginary numbers if you feed the random
%               gaussian to them (that is the energy) with the negative values
%               still intact, to compensate for this we feed them the absolute
%               value of the energy and then choose a random direction, plus or
%               minus, for that offset to move in; also, the heads that would
%               previously have ghosted past the actin filament no longer do so
% CDW(20080131)-The randn now hits the distance moved instead of the energy
%               causing the movement, this seems wrong and is not as nice
%               as the mechanism before, but it produces a monomial
%               distribution of locations instead of the weird binomial in
%               each dimension that we were getting
% CDW(20080521)-The randn now hits a force. The springs move until this force is
%               countered by the force they exert by their displacement 


%% Code

%Uncomment this next line when a major new ver is saved
% if ~strcmp(lastwarn, ['Running an old version of ' mfilename]) %only warn once
%   warning('DangerDave:OldVersion',['Running an old version of ' mfilename]) 
% end

%find all the unbound heads (their linear indices)
UnbHds = find(Mf.bst == 0);

%% Get the forces we will use
ForceMag = randn(1, length(UnbHds)).*Sc.dfv;
ForceDir = rand(1,length(UnbHds)).*2*pi;
ForcePerp = ForceMag .* sin(Mf.ths(UnbHds) - ForceDir);
ForcePara = ForceMag .* cos(Mf.ths(UnbHds) - ForceDir);

%% Find new location offsets in pol and cartesian coords
Pol(1,:) = (ForcePara / Mf.rk(UnbHds)) + Mf.rs(UnbHds);
Pol(2,:) = ((ForcePerp .* Pol(1,:)) ./ Mf.thk(UnbHds)) + Mf.ths(UnbHds);
[Cart(1,:), Cart(2,:)] = pol2cart(Pol(2,:), Pol(1,:));

%% Update the Mf head locations
Mf.hloc(1,UnbHds) = Cart(1,:) + Mf.loc(1,UnbHds); %update the head location 
Mf.hloc(2,UnbHds) = Cart(2,:) + Mf.loc(2,UnbHds); %from the calculated offsets
Mf.hloc(2,Mf.hloc(2,UnbHds)>Sc.sep) = Sc.sep; %stop heads from ghosting past the thin filament



%% Previous code
%CDW20080521 - Code has been refactored so many times that it is now unreadable
%              for the most part. As such it has been archived below and the
%              quest continues above
% %and create some local variables
% T = 288;                %the temperature (in K) that this runs at
% K = 1.381 * 10^-23;     %Boltzman const (in J/K)
% 
% %find all the unbound heads (their linear indices)
% UnbHds = find(Mf.bst == 0);
% 
% %%create a matrix that will contain our random energies and directions thereof in
% %%the manner: [r,...;th,...]
% DiffusionPos = zeros(2,size(UnbHds,2));
% 
% % %Determine the radial offset from the myosin node that each of our unbound heads
% % % should possess
% % RandNum = rand(1,length(UnbHds))>.5;  %to make the gaussing profile of the energy go both ways (without throwing 
% % RandDir(RandNum==0)=-1;               %imaginary components as it would do in the abscence of an abs()  )
% % RandDir(RandNum==1)=1;
% % DiffusionPos(1,:) = sqrt(2                   *K*T.*abs(randn(size(UnbHds)))./Mf.rk(UnbHds)).* sqrt(10^21).*RandDir + Mf.rs(UnbHds); %the sqrt(10^21) is a conversion factor to make the term before it convert over to nm
% DiffusionPos(1,:) = .7*randn(size(UnbHds)).*sqrt(2                        *K*T./Mf.rk(UnbHds)).* sqrt(10^21) + Mf.rs(UnbHds); %the sqrt(10^21) is a conversion factor to make the term before it convert over to nm
% 
% 
% % %We (unfortunetly) need the current value of r to calculate the proper Th
% % [CurrTh, CurrR] = cart2pol(Mf.hloc(1, UnboundHeads)-Mf.loc(1, UnboundHeads),Mf.hloc(2, UnboundHeads)-Mf.loc(2, UnboundHeads));
% % CurrR = DiffusionPos(1,:) + CurrR;
% 
% % %Determine the angle that each of our unbound myosin heads should be located at,
% % % relative to the myosin node and the x axis
% % RandNum = rand(1,length(UnbHds))>.5;  %to make the gaussing profile of the energy go both ways (without throwing 
% % RandDir(RandNum==0)=-1;               %imaginary components as it would do in the abscence of an abs()  )
% % RandDir(RandNum==1)=1;
% % DiffusionPos(2,:) = sqrt(2.*DiffusionPos(1,:)*K*T.*abs(randn(size(UnbHds)))./Mf.thk(UnbHds)).*sqrt(10^21).*RandDir + Mf.ths(UnbHds);
% DiffusionPos(2,:) = .7*randn(size(UnbHds)).*sqrt(2.*abs(DiffusionPos(1,:)).*K.*T./Mf.thk(UnbHds)).*sqrt(10^21) + Mf.ths(UnbHds);
% 
% 
% %Copy our nice new values back to the Mf.hloc matrix we will be returning in the
% % Mf structure
% [Mf.hloc(1,UnbHds), Mf.hloc(2,UnbHds)] = pol2cart(DiffusionPos(2,:), DiffusionPos(1,:));
% Mf.hloc(1:2,UnbHds) = Mf.loc(1:2,UnbHds)+Mf.hloc(1:2,UnbHds);
% 
% %And any of the fellows that got too big for their britches and ghosted past the
% %thin filament are snapped back to the thin filament y location, 
% %this occurs by looking at all the unbound nodes, multiplying those that are not
% %over the line by one and those that are over the line are set to one and 
% %multiplied by the line
% Mf.hloc(2,UnbHds) = (Mf.hloc(2,UnbHds)<Sc.sep).*Mf.hloc(2,UnbHds)+(Mf.hloc(2,UnbHds)>Sc.sep).*Sc.sep;
% 
% % %flag those entries as the ones getting a force exerted on them or an angle
% % %selected for them
% % DiffusionF(1, UnboundHeads) = 1;
% % DiffusionF(2:3,UnboundHeads) = 2*pi;
% % 
% % %generate the random forces and random angles
% % DiffusionF(1,:) = DiffusionF(1,:) .* normrnd(ForceMean, ForceVar, [1,size(DiffusionF,2)]);
% % DiffusionF(2:3,:) = DiffusionF(2:3,:) .* rand([2,size(DiffusionF,2)]);
% % 
% % %calculate offsets as a result of the forces
% % 
% % Offsets = motorF2Pos_Ver1(DiffusionF,Mf);
% % 
% % %update the thick filament
% % Mf.hloc = Mf.hloc + [Offsets];