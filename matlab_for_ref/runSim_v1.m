function [Log, InitialConditions]=runSim_v1(Mf, Af, Sc, Runs, Iters, Path)

% CDW 20080117 - This is it, the lovely fellow that is going to take in some
%                inital conditions and information about how many times to look
%                over them and how long to do it each time, and return to you a
%                log of all that occured and the initial conditions again, I
%                haven't yet decided if it should save the output to disk or let
%                the parent script do that, but I am leaning towards letting the
%                parent do that
% CDW 20080521 - The parent does that

%% Old verision protection
%Uncomment this next line when a major new ver is saved
% if ~strcmp(lastwarn, ['Running an old version of ' mfilename]) %only warn once
%   warning('DangerDave:OldVersion',['Running an old version of ' mfilename]) 
% end


%% Initialize bit and pieces

%initialize our random number generators, very important to make sure that all
%our trials look different
randn('state', sum(100*clock));
rand('state' , sum(100*clock));

%initalize the function that holds all the exit flag info
exitflag = NaN(2,Iters);

%initialize the output logging structure
Log.atp   = zeros(Iters, length(Mf.atp));
Log.bst   = zeros(Iters, length(Mf.bst));
Log.axial = zeros(Iters, 2);
Log.perp  = zeros(Iters, length(Mf.loc));

%initialize the structure of initial conditons we will return, this is helpful
%as we may be calling this function a good many times with different inputs
InitialConditions.Mf    = Mf;
InitialConditions.Af    = Af;
InitialConditions.Sc    = Sc;
InitialConditions.Runs  = Runs;
InitialConditions.Iters = Iters;
InitialConditions.Path  = Path;

%initalize a structure of known locations
KnownLoc = struct('bst', Mf.bst, 'aloc', Af.loc, 'mloc', Mf.loc, 'hloc', Mf.hloc);


%% Perform the runs

for CurrRun = 1:Runs
    
    %re-initilaize our filaments for each run
    Mf = InitialConditions.Mf;
    Af = InitialConditions.Af;
    Sc = InitialConditions.Sc;
    
    for CurrIter = 1:Iters
        %remember our current head state
        HeadState = Mf.bst;

        %perturb the motor head locations
        [Mf] = motorPerturb_Ver3(Mf, Sc);

        %transition motors from one bound state to another (maybe)
        [Mf,Af] = motorBind_Ver2(Mf, Af, Sc);

        %balance the forces from any state transitions
        %if something has changed,
        if ~isempty(find(Mf.bst-HeadState, 1)),
            %and we haven't seen it before, optimize
            if isempty(find(arrayfun(@(x)(isequal(x.bst,Mf.bst)),KnownLoc), 1))
                %by sending the filaments to a balance function that calls fsolve
                [Mf, Af, exitflag(1,CurrIter)] = balFils_v4(Mf,Af,Sc,'gn');
                %if that didn't get us close to balanced, call another type of
                %optimization function that get's us to a minima in a different way
                if exitflag(1,CurrIter) ~= 1,
                    %specifically, using a dogleg method instead of a gauss-newton
                    [Mf, Af, exitflag(2,CurrIter)] = balFils_v4(Mf,Af,Sc,'dogleg');
                end
                %and update the known locations
                KnownLoc(end+1).bst = Mf.bst; %this guy adds a new entry into the KnownLoc strcture array as well
                KnownLoc(end).mloc = Mf.loc; 
                KnownLoc(end).hloc = Mf.hloc;
                KnownLoc(end).aloc = Af.loc;
            else %when we have seen this config before
                KnownCase = find(arrayfun(@(x)(isequal(x.bst,Mf.bst)),KnownLoc)); %figure out which config it is
                Mf.bst = KnownLoc(KnownCase).bst; %and copy the solution back to our fil data structures
                Mf.loc = KnownLoc(KnownCase).mloc; 
                Mf.hloc = KnownLoc(KnownCase).hloc;
                Af.loc = KnownLoc(KnownCase).aloc;
            end
        end

        %log the output
        Log = logOutput_Ver1(Mf, Af, Log, CurrRun, CurrIter) ;
        %reset the atp utilization record so we don't see this turn's atp useage again
        Mf.atp = zeros(1,length(Mf.atp));

    end

%     %tell us where we are in the total number of runs we have left
%     if mod(run,1)==0
%        	disp(['On ' num2str(run) ' of ' num2str(runs)])
%     end    
end

%save(['sep' num2str(Sc.sep) '.mat'], 'Log')

%disp(['Took ' num2str(t,3) ' seconds to run the whole file'])
%disp(['Node(s) ' num2str(find(Mf.bst)) ' were bound to node(s) ' num2str(Mf.lnk(find(Mf.bst)))])


% if isunix
%  !/usr/bin/mail cdave@u.washington.edu -s 'Cluster job complete' < /home/cluster/thedrick/message.txt
% end

% catch
%     %this throws an email to me if there is an error while the file is running
%     if isunix
%         ! echo 'An error has occured, the cluster needs looking to.' | /usr/bin/mail cdave@u.washington.edu -s 'Cluster Error' 
%     end
% end


