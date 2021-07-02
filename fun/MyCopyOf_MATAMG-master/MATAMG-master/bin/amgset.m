function options = amgset(varargin)
%AMGSET Create/alter amg OPTIONS structure.
%  OPTIONS = AMGSET('PARAM1',VALUE1,'PARAM2',VALUE2,...) create an
%  amg option structure OPTIONS in which the named parameters have
%  the specified values. Any unspecified parameters are set to [](parameters
%  with value [] indicate to use the default value for the parameter when
%  OPTIONS is passed to the amg function). It is sufficient to type
%  only the leading characters that uniquely identify the parameter. Case is
%  ignored for parameter names.
%  NOTE: For values that are strings, the complete string is required.
%
%  OPTIONS = AMGSET(OLDOPTS,'PARAM1',VALUE1) creates a copy of OLDOPTS
%  with the named parameters altered with the specified values.
%
%  OPTIONS = AMGSET(OLDOPTS,NEWOPTS) combines an existing options structure
%  OLDOPTS with a new options structure NEWOPTS
%  with non-empty values overwirte the corresponding old parameters in
%  OLDPLOTS.
%
%  AMGSET with no input arguments and no output arguments displays all
%  parameter names and their possible values, with defaults shown in {}
%  when the default is the same for all functions that use that parameter.
%
%  OPTIONS = AMGSET(with no input arguments) creates an options structure
%  OPTIONS where all the fields are set to [].
%
%AMGSET PARAMETERS for MATLAB
%
%  Examples:
%    
%
%  See also amgget.

%   Copyright 2011 Minho Park.
%   
% Changelog
%   16-05-2011 Added aggressive coarsening and interpolation
%   v.0.1.2 - Added aggressive coarsening and interpolation (A(2,1) A(2,2))

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('%15s\n','[AMG]');
    fprintf('%15s : %s\n','Coarsest','[ positive scalar | {2} ]');
    fprintf('%15s : %s\n','MaxCycle','[ positive scalar | {20} ]');
    fprintf('%15s : %s\n','TolAMG','[ positive scalar | {1e-10} ]');
    fprintf('%15s : %s\n','Iter1','[ positive scalar | {1} ]');
    fprintf('%15s : %s\n','Iter2','[ positive scalar | {1} ]');
    fprintf('%15s : %s\n','MuCycle','[ positive scalar | {1} ]');
    
    fprintf('%15s\n','[INTERPOLATION]');
    fprintf('%15s : %s\n','IntpType','[ {amg} | aamg | rbamg | lramg | off ]');
    fprintf('%15s : %s\n','NumVec','[ positive scalar | 5 ]');
    fprintf('%15s : %s\n','NumRel','[ positive scalar | 5 ]');
    
    fprintf('%15s\n','[RELAXATION]');
    fprintf('%15s : %s\n','RelType','[ {gs} | jacobi ]');
    fprintf('%15s : %s\n','RelPara','[ positive scalar | {1.0} ]');
    
    fprintf('%15s\n','[COARSENING]');
    fprintf('%15s : %s\n','CsnType','[ {amg} | amg(abs) | cr | load | a21 | a22 ]');
    fprintf('%15s : %s\n','Theta','[ positive scalar | {0.25} ]');
    fprintf('%15s : %s\n','Alpha','[ positive scalar | {0.7} ]');
    fprintf('%15s : %s\n','Nu','[ positive scalar | {5} ]');
    
    fprintf('%15s\n','[   FILE   ]');
    fprintf('%15s : %s\n','Log','[ on | {off} ]');
    fprintf('%15s : %s\n','LogFile','[ {out} ]');
    fprintf('%15s : %s\n','LogNum','[ {on} | off ]');
    fprintf('%15s : %s\n','SaveCsn','[ on | {off} ]');
    fprintf('%15s : %s\n','SaveIntp','[ on | {off} ]');
    
    fprintf('%15s\n','[ PRECONDITIONER ]');
    fprintf('%15s : %s\n','PreCond','[ {off} | pcg ]');
    
    fprintf('%15s\n','[ MISC ]');
    fprintf('%15s : %s\n','PrintOnScreen','[ {on} | off ]');
    return;
end

% Create a cell array of all the field names
allfields = {'Coarsest';'MaxCycle' ;'TolAMG';'Iter1';'Iter2';...
    'MuCycle';'IntpType';'NumVec';'NumRel';'RelType';'RelPara';...
    'CsnType';'Theta';'Alpha';'Nu';'Log';'LogFile';'LogNum';'SaveCsn';...
    'SaveIntp';'PreCond';'PrintOnScreen';'Numbering'};

% Create a struct of all the field with all values set to []
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
%[]'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assing it

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    
    if ischar(arg)  % arg is an option name
        break;
    end
    if ~isempty(arg)    % [] is a valid options argument
        if ~isa(arg,'struct')
            error('MATLAB:amgset:NoParamNameOrStruct',...
                ['Expected argument %d to be a string parameter name '...
                'or an options structure \ncreated with AMGSET.'],i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% finite to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('MATLAB:amgset:ArgNameValueMismatch',...
        'Arguments must occur in name-value paris.');
end
    
expectval = 0;      % start expecting a name, not a valu

while i <= numberargs
    
    arg = varargin{i};
    if ~expectval
        if ~ischar(arg)
            error('MATLAB:amgset:InvalidParamName',...
                'Expected argument %d to be a string parameter name.',i);
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        
        if isempty(j)       % if no matches
            % Error out -compose internationalization-friendly message with
            % hyperlinks
            stringWithLink = formatStringWithHyperlinks(sprintf('Link to reference page'),'doc amgset');
            error('MATLAB:amgset:InvalidParamName',...
                ['Unrecognized parameter name ''%s''. Please see the optimset' ...
                   ' reference page in the documentation for a list of acceptable' ...
                   ' option parameters. %s'],arg,stringWithLink);
        elseif length(j) > 1        % if more than one match
            % Check for any exact matches (in case any names are subnsets
            % of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambigous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ',' Names{k,:}];
                end
                msg = [msg,'.'];
                error('MATLAB:amgset:AmbiguousParamName',msg);
            end
        end
        expectval = 1;  % we expect a value next
        
    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        
        checkfield(names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i= i + 1;
end

if expectval
    error('MATLAB:amgset:NoValueForParam',...
        'Expected value for parameter ''%s''.',arg);
end

function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V,OPTIMTBX) checks the contents of the specified
%   value V to be valid for the field 'field'. OPTIMTBX indicates if 
%   the Optimization Toolbox is on the path.
%

% empty matrix is always valid
if isempty(value)
    return
end

field = lower(field);
% See if it is one of the valid MATLAB fields.  It may be both an Optim
% and MATLAB field, e.g. MaxFunEvals, in which case the MATLAB valid
% test may fail and the Optim one may pass.
validfield = true;

switch field
    case {'coarsest'} % integer scalar >= 2
        [validvalue, errmsg, errid] = IntGTE(field,value,2);
    case {'maxcycle','mucycle','numvec','numrel','nu'} % integer scalar >= 1
        [validvalue, errmsg, errid] = IntGTE(field,value,1);
    case {'iter1','iter2'}  % integer scalar >= 0
        [validvalue, errmsg, errid] = IntGTE(field,value,0);
    case {'tolamg','relpara','theta','alpha'} % positive real scalar
        [validvalue, errmsg, errid] = posReal(field,value);
    case {'savecsn','saveintep','printonscreen','log','lognum'} % off,on
        [validvalue, errmsg, errid] = onOffType(field,value);
    case {'csntype'} % several character string for coarsening option
        [validvalue, errmsg, errid] = csnType(field,value);
    case {'reltype'}    % several character string for relaxation option
        [validvalue, errmsg, errid] = relType(field,value);
    case {'intptype'}    % several character string for relaxation option
        [validvalue, errmsg, errid] = intpType(field,value);
    case {'precond'}    % several character string for relaxation option
        [validvalue, errmsg, errid] = preCondType(field,value);
    case {'logfile'}
        if ~ischar(value)
            errmsg = 'MATLAB:funfun:amgset:LogFile:NotAString';
            errid = 'Invalid value for OPTIONS parameter LogFile: must be a string';
            validvalue = false;
        else
            errid = '';
            errmsg = '';
            validvalue = true;
        end            
    otherwise
        validfield = false;
        validvalue = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
        errid = 'MATLAB:optimset:checkfield:InvalidParamName';
end

if validvalue 
    return;
elseif validfield  
    % Throw the MATLAB invalid value error
    error(errid, errmsg);
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = IntGTE(field,value,lowerbnd)
% Any nonnegative real scalar or sometimes a special string

valid =  isreal(value) && isscalar(value) && (value >= lowerbnd) && value == floor(value);

if ~valid
    if ischar(value)
        errid = 'MATLAB:funfun:amgset:IntGTl:IntLTNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a integer scalar >=  %d (not a string).',field,lowerbnd);
    else
        errid = 'MATLAB:funfun:amgset:IntGTl:IntLTNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a integer scalar >=  %d.',field,lowerbnd);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = posReal(field,value)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value > 0) ;

if ~valid
    if ischar(value)
        errid = 'MATLAB:funfun:amgset:posReal:nonnegativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar (not a string).',field);
    else
        errid = 'MATLAB:funfun:amgset:posReal:nonnegativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = csnType(field,value)
% One of these strings: amg, cr, load
valid =  ischar(value) && any(strcmp(value, ...
    {'amg';'amg(abs)';'cr';'load';'a21';'a22'}));
if ~valid
    errid = 'MATLAB:funfun:amgset:csnType:notACsnType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s: must be ''amg'',''amg(abs)'',''cr'',''a21'', ''a22'' or ''load''.'],field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = relType(field,value)
% One of these strings: gs, jacobi
valid =  ischar(value) && any(strcmp(value, ...
    {'gs';'jacobi'}));
if ~valid
    errid = 'MATLAB:funfun:amgset:relType:notARelType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s: must be ''gs'' or ''jacobi''.'],field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = intpType(field,value)
% One of these strings: amg, aamg, rbamg, off
valid =  ischar(value) && any(strcmp(value, ...
    {'amg';'aamg';'rbamg';'lramg';'off'}));
if ~valid
    errid = 'MATLAB:funfun:amgset:intpType:notAIntpType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s: must be ''amg'', ''aamg'', ''rbamg'',''lramg'' or ''off''.'],field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = preCondType(field,value)
% One of these strings: off, pcg
valid =  ischar(value) && any(strcmp(value, ...
    {'off';'pcg'}));
if ~valid
    errid = 'MATLAB:funfun:amgset:preCondType:notAPreCondType';
    errmsg = sprintf(['Invalid value for OPTIONS parameter %s: must be ''off'' or ''pcg''.'],field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
    errid = 'MATLAB:funfun:optimset:onOffType:notOnOffType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function formattedString = formatStringWithHyperlinks(textToHyperlink,commandToRun)
% Check if user is running MATLAB desktop. In this case wrap
% textToHyperlink with HTML tags so that when user clicks on
% textToHyperlink, commandToRun gets executed. If not running
% MATLAB desktop, leave textToHyperlink unchanged.

if feature('hotlinks') && ~isdeployed
    % If using MATLAB desktop and not deployed, use hyperlinks
    formattedString = sprintf('<a href="matlab: %s ">%s</a>.',commandToRun,textToHyperlink);
else
    % Use plain string
    formattedString = sprintf('');
end


             