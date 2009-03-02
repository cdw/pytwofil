function [jm, jh] = startClusterJob_v3(Mf, Af, Sc, Runs, Iters, Param)

%Version 1:
%This function exists to start a cluster job, it takes as it's input the initial
%filaments the cluster job should start with and then parses out the runs to the
%workers such that they are roughly evenly distributed and each worker has to do
%two (or one?) tasks before the job is complete. It returns the job handle (jh) 
%and the job manager (jm) and, finally, lets the calling script/function do the 
%remaining work. - CDW 20080117
%Version 2:
%This now uses the 'FinishedFcn' property of a job to call whenJobIsFinished_v1 
%at the end of each job, removing the need for this to be done by some
%other orgainizing function - CDW 20080119
%Version 3:
% This now takes in Param, which has Param.name and Param.value which let
% us name the tasks in a more helpful fashion

%% Set variables we are hard coding
%this is the location that we are using to set variables that in theory could be
%passed into this function but which I am instead hard-coding into the function
%in order to simplify things
% PathToCode = {'/exports/cluster/dave/matlab/bps2008/code'}; %where the tasks can find the code they'll need
% PathToData = ['/exports/cluster/dave/matlab/bps2008/data/' datestr(date, 'yyyymmdd')];  %where we will want to store returned data
[PathToCode, PathToData] = getPaths_v1();

%% Old verision protection
%Uncomment this next line when a new ver is saved
% if ~strcmp(lastwarn, ['Running an old version of ' mfilename]) %only warn once
%   warning('DangerDave:OldVersion',['Running an old version of ' mfilename]) 
% end


%% Find the job manager 
% this can be in the local or the cluster case, although really it will only
% work in the cluster case for now, as a function of relying on the
% NumerOfIdleWorkers property of a real proper job manager, this can be fixed by
% making the sections below this one diffrent in the pc and unix cases - CDW
% 20080117

if ispc==1 %if on my laptop
    jm = findResource('scheduler','type','local');
else %if on the cluster
    jm = findResource('jobmanager','Name','jm_bert','LookupURL','192.168.99.1');
end


%% Set up needed inputs to the job
% that is, figure out the task number and what have you

IdleWorkers     = get(jm, 'NumberOfIdleWorkers');
BusyWorkers     = get(jm, 'NumberOfBusyWorkers');
TotalWorkers    = BusyWorkers + IdleWorkers;
TaskWorkers     = 15; % This is here to make it easy to switch to only 
                                % some workers if needed now we know the numbers 
                                % of all the workers, set up the number per task
% Our goal, again, is to have each worker complete one task
RunsPerTask     = repmat(floor(Runs./TaskWorkers),1,TaskWorkers);
RunsPerTask(1)  = RunsPerTask(1) + mod(Runs,TaskWorkers);
% Now we need to make the cell array of input cells
InputCells      = cell(1,TaskWorkers);
for i = 1:TaskWorkers,
    InputCells{i} = {Mf, Af, Sc, RunsPerTask(i), Iters, PathToData};
end


%% Start the job

%Create it
jh = createJob(jm);
%Set some properties
set(jh, 'PathDependencies', {PathToCode}); %tell it where to find the code it will need
set(jh, 'Name', ['dave_' Param.name '_' sprintf('%07.2f',Param.value)]);%datestr(now, 'yymmddtHH:MM:SS.FFF')]); %Changed name to use sep as that is really the identifing characteristic
set(jh, 'FinishedFcn', @whenJobIsFinished_v1); %give a try to letting jobs process themselves when done
%Populate it
for TaskNum = 1:length(InputCells)
  createTask(jh, @runSim_v1, 2, InputCells{TaskNum});  %create tasks that reside within the job
end

%Don't know if the next bit will work, am just gonna have script wait until job
%is done for the moment
% %Make it send an email and clean up after itself when done
% set(jh, 'FinishedFcn', @whenJobIsFinished_v1);

%Get it running
submit(jh);









