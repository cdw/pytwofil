function[Mf, Af, Sc]=initFils(MEnd, AEnd)

%[Mf, Af]=initFils(MZeroLoc, AZeroLoc)
%
% this initializes one thick and one thing filament and returns two 
% structures containing the filaments and their associated properties
%
% Mf - structure of initialized thick filament
% Af - structure of initialized thin filament
% Sc - structure of sarcomere properties that don't fit into Mf or Af
% MEnd - 3 vector of the end attachment point of the thick filament
% AEnd - 3 vector of the end attachment point of the thin filament
%
%
% Output structures
% Myosin Filament
% Mf.loc	(x1,...;y1,...;z1,...) 	-cartesian coordinates of the thick fil node
% 	.hloc	" "                     -cartesian coordinates of the motor head
% 	.bst	(0,2,1,...)             -binding state of the node's motor
% 	.lnk	(0,12,14,...)           -actin node motor is bound to (or 0)
% 	.rs		(2,4,2,...)             -rest length of motor linear spring
% 	.rk		" "                     -spring constant of motor linear spring
% 	.rv		(rs1,rs2;rk1,rk2)		-values of rs, rk for different bound states
% 	.ths	(2,4,2,...)             -rest angle \theta of motor
% 	.thk	" "                     -rot. spring constant in \theta dir
% 	.thv	(ths1,ths2;thk1,thk2)	-values of ths, thk for different bound	states
% 	.phs	(2,4,2,...)             -rest angle \phi of motor
% 	.phk	" "                     -rot. spring constant in \phi dir
% 	.phv	(phs1,phs2;phk1,phk2)	-values of phs, phk for different bound states
% 	.tis	??						-rest length of connective titin
% 	.tik	??						-spring constant of connective titin
% 	.ms		43						-rest length between thick fil nodes
% 	.mk		2020					-spring const of the thick fil
% 	.uda	40						-rest length of undecorated area
%   .mln    (0,0,0)                 -m line attachment point
%   .zln    (1000,3,0;1000,-3,0)    -z line attachment points (of titin)
% Actin Filament
% Af.loc	(x1,...;y1,...;z1,...)	-cartesian coordinates of the thin fil node
% 	.bst	" "                     -binding state of node
% 	.lnk	(0,2,1,...)             -motor the node is linked to (or 0)
% 	.as     37						-rest length between nodes 
% 	.ak		1745 					-spring constant of the thin fil 
%   .zln    (1000,11.3,0)           -z line attachment point
% Sarcomere
% Sc.len 	1400                    -the total length of the half sarcomere
%   .sep 	11.3                    -the distance between thick and thin fils
%   .dfm 	0                       -the mean force of diffusion
%   .dfv    2                       -the variance of the diffusional force
%   .bd     .55                     -the dist at which a motor and actin bind

%set up variables
Mk = 2020;      %spring const between thick fil sites
Ms = 43;        %rest length between thick fil sites
Mn = 20;        %number of thick fil sites
Ak = 1743;      %spring const between actin sites
As = 37.3;      %rest length between thin fil sites
An = 30;        %number of thin fil sites
%Tik = 2020;     %the spring constant of titin
Uda = 40;       %how far last thin site is from first thick site at rest
Ths = pi/3;     %rest angle theta
Phs = pi/2;     %rest angle phi
Rs = 10.5;      %crossbridge restlength
Rk = 5;         %spring const of cross bridges
Thk = 200;       %spring const of torsional spring theta
Phk = 80;       %spring const of torsional spring phi
%Sep = 11.3;     %distance between thick and thin filaments
%TitinOff = 2;   %distance titin is offset from being in line with the thick fil
ForceMean = 0;  %mean of forces exerted on myosin heads
ForceVar = 10.0; %variance of forces exerted on myosin heads
BindDist= .55; %.55   %the distance from a binding site at which you are considered bound


% initialize the location of all the sites where the myosin tails are bound
Mf.loc  = [MEnd(1)+Uda:Ms:(MEnd(1)+Uda+(Mn-1)*Ms); ...
                repmat(MEnd(2), 1, Mn); ...
                repmat(MEnd(3), 1, Mn)];
Mf.mln  = MEnd;
%No Titin for now - CDW 20080114
% Mf.zln  = [AEnd(1),MEnd(2)+TitinOff,MEnd(3); ...
%            AEnd(1),MEnd(2)-TitinOff,MEnd(3)];

% initialize the values for the myosin heads themselves
Mf.hloc = Mf.loc+... 
                repmat([Rs*cos(Ths);Rs*sin(Ths);0], 1, Mn);
Mf.bst  = repmat(0, 1, Mn);
Mf.lnk  = repmat(0, 1, Mn);
Mf.atp  = repmat(0, 1, Mn);

% initialize all the various constants and variables we will refer to 
Mf.rs  = repmat(Rs,  1, Mn);
Mf.rk  = repmat(Rk,  1, Mn);
Mf.ths = repmat(Ths, 1, Mn);
Mf.thk = repmat(Thk, 1, Mn);
Mf.phs = repmat(Phs, 1, Mn);
Mf.phk = repmat(Phk, 1, Mn);
%Mf.tis = norm([(MEnd(1)+Uda+Ms*(Mn-1)),MEnd(2)]-[Mf.zln(1,1),Mf.zln(1,2)]); %dist from last thick node to a titin bind site
%Mf.tik = Tik;
Mf.ms  = Ms ;
Mf.mk  = Mk ; 
Mf.uda = Uda;

Mf.rv  = [Rs,Rs,0.9*Rs    ;Rk,Rk,Rk]; %[unbound, loosely bound, tightly bound]
Mf.thv = [Ths,Ths,1.2*Ths ;Thk,Thk,Thk];
Mf.phv = [Phs,Phs,Phs     ;Phk,Phk,Phk];


Af.loc =  [(AEnd(1)-An*As):As:AEnd(1)-As; ...
                repmat(AEnd(2), 1, An); ...
                repmat(AEnd(3), 1, An)];
Af.zln = AEnd;
Af.bst = repmat(0, 1, An);
Af.lnk = repmat(0, 1, An);
Af.as  = As;
Af.ak  = Ak;
	
Sc.len = AEnd(1) - MEnd(1);
Sc.sep = AEnd(2) - MEnd(2);
Sc.dfm = ForceMean;
Sc.dfv = ForceVar;
Sc.bd  = BindDist;
