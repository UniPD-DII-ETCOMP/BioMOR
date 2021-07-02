function o = amgget(options,name,default)
%AMGGET Get amg OPTIONS parameters
%  VAL = AMGGET(OPTIONS,'NAME') extracts the value of the named parameter
%  from amg options structure OPTIONS, returning an empty matrix if 
%  the parameter value is not specified in OPTIONS. it is sufficient to
%  type only the leading characters that uniquly identify the 
%  parameter. Case is ignored for parameter names. [] is a valid OPTIONS
%  argument.
%
%  VAL = AMGGET(OPTIONS,'NAME',DEFAULT) extract the named parameter as
%  ab-ove, but returns DEFAUT if the named parameter is not specified (is [])
%  in OPTIONS. For example
%
%     val = amgget(options,'CsnType','cr');
%
%  returns val = 'cr' if the CsnType parameter is not specified in options.
% 
%  See also amgset.

% CHANGE
%   13-05-2011 Added properties PRINTONSCREEN AND NUMBERING


if nargin < 2
    error('MATLAB:amgget:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('MATLAB:amgget:Arg1NotStruct',...
        'First argument must be an options structure created with AMGSET.');
end

if isempty(options)
    o = default;
    return;
end

% Create a cell array of all the field names
allfields = {'Coarsest';'MaxCycle';'TolAMG';'Iter1';'Iter2';...
    'MuCycle';'IntpType';'NumVec';'NumRel';'RelType';'RelPara';...
    'CsnType';'Theta';'Alpha';'Nu';'Log';'LogFile';'LogNum';'SaveCsn';'SaveIntp';...
    'PreCond';'PrintOnScreen';'Numbering'};

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));

if isempty(j)               % if no matches
    error('MATLAB:amgget:InvalidPropName',...
        ['Unrecognized option name ''%s''.  ' ...
        'See AMGSET for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous option name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = [msg, '.)'];
        error('MATLAB:optimget:AmbiguousPropName', msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    o = options.(Names{j,:});
    if isempty(o)
        o = default;
    end
else
    o = default;
end
