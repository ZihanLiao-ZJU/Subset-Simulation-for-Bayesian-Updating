function core=setup_semcs(core)
% Set SEMCS_ROOT for a development checkout. Release dependency is SHA pinned.
% Select the core before creating any sampler/model objects in a fresh session.
% Loaded classes are never automatically cleared to switch core versions.
root=fileparts(mfilename('fullpath'));
if nargin<1||isempty(core),core=getenv('SEMCS_ROOT');end
if isempty(core)
 if exist(fullfile(root,'+semcs','Engine.m'),'file'),core=root;
 else
  pin=jsondecode(fileread(fullfile(root,'semcs_dependency.json')));
  assert(~strcmp(pin.commit,'PENDING_CANONICAL_RELEASE'),'semcs:Unpinned','Release dependency commit has not been set.');
  cache=fullfile(root,'.semcs-deps',pin.commit);core=fullfile(cache,['Sequential-multiple-importance-sampling-for-Reliability-Estimation-' pin.commit]);
  if ~exist(fullfile(core,'+semcs','Engine.m'),'file')
   if ~exist(cache,'dir'),mkdir(cache);end
   archive=fullfile(cache,'core.zip');url=['https://github.com/ZihanLiao-ZJU/Sequential-multiple-importance-sampling-for-Reliability-Estimation/archive/' pin.commit '.zip'];
   websave(archive,url);unzip(archive,cache);delete(archive);
  end
 end
end
assert(exist(fullfile(core,'+semcs','Engine.m'),'file')~=0,'semcs:Dependency','SEMCS_ROOT must contain +semcs/Engine.m');
[valid,attributes]=fileattrib(core);
assert(valid&&isscalar(attributes)&&attributes.directory,'semcs:Dependency','Core root must be one existing directory.');
core=attributes.Name;
addpath(core,'-begin');
resolved=which('semcs.Engine');expected=fullfile(core,'+semcs','Engine.m');
resolved=strrep(resolved,'\','/');expected=strrep(expected,'\','/');
if ispc,same=strcmpi(resolved,expected);else,same=strcmp(resolved,expected);end
assert(same,'semcs:CoreShadowed', ...
    'Selected core does not resolve semcs.Engine. Start a fresh MATLAB session before selecting a different core.');
end
