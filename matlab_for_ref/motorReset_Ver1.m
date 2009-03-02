function[Mf]=motorReset_Ver1(Mf)

% [Mf]=resetUnboundMHeads_Ver1(Mf)
%
% this resets all the unbound myosin heads to their rest position
%
% Mf            - structure containing all info about a filament
% Mf.thv(1,1)   - the rest angle from the thick filament to the myosin body
% Mf.rv(1,1)    - the rest length of the cross bridge
% [Mf]          - the updated filament

%% General Documentation
% CDW(20070731)-This is one of the retrofits from the rotation code


%% Code

%Uncomment this next line when a new ver is saved
%warning(['Running an old version of ' mfilename]) 

Rs = Mf.rv(1,1);
Ths = Mf.thv(1,1);
Phs = Mf.phv(1,1);

% capture the rest angle and length as an offset in cartesian coordinates
RestOffset = [Rs*sin(Phs)*cos(Ths); Rs*sin(Phs)*sin(Ths); Rs*cos(Phs)];

% find out which heads are unbound, they are the ones we will work with
UHeads = find(Mf.bst == 0);

% set all unbound heads to be the corresponding thick fil loc plus offset
Mf.hloc(:, UHeads) = ...
    Mf.loc(:, UHeads) + repmat(RestOffset,1, length(UHeads));

