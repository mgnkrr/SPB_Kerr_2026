function S = add_g_fields(~, varargin)
% Add Gmin, optional advection correction, and per-model GΔ into S

p = inputParser;
addParameter(p,'datafile',fullfile('datasets','gmin_data.mat'),@(s)ischar(s)||isstring(s));
addParameter(p,'advection_file','datasets/longitudinal_advection_maps.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'advection_var','Delta_q_adv',@(s)ischar(s)||isstring(s));
addParameter(p,'advection_units','mWm2',@(s)any(strcmpi(s,{'mWm2','Wm2'})));
addParameter(p,'keep_gdiff_noadv',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'cast_like_H',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'save',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'outfile',fullfile('datasets','gmin_data.mat'),@(s)ischar(s)||isstring(s));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
opt = p.Results;

logf = @(varargin) ( opt.verbose && fprintf('[add_g_fields] %s\n', sprintf(varargin{:})) );

% ----------- load S -----------
datafile = char(opt.datafile);
assert(isfile(datafile), 'Data file not found: %s', datafile);
L = load(datafile,'S'); S = L.S; clear L;

% ----------- sanity -----------
mustHave = {'Xgrid','Ygrid','H','Ts','Mb','names','models'};
for k = 1:numel(mustHave)
    assert(isfield(S, mustHave{k}), 'S.%s missing in %s', mustHave{k}, datafile);
end

% Ensure model fieldnames are valid & unique 
validNames = cellfun(@matlab.lang.makeValidName, S.names, 'uni', 0);
validNames = matlab.lang.makeUniqueStrings(validNames);

% Verify models actually have these fields
modelFields = fieldnames(S.models);

% Attempt to align by valid name, else by case-insensitive match
for i = 1:numel(validNames)
    vn = validNames{i};
    if ~isfield(S.models, vn)
        % try to locate a close match
        hit = find(strcmpi(modelFields, vn), 1);
        if ~isempty(hit)
            % rename model field to vn
            old = modelFields{hit};
            S.models.(vn) = S.models.(old);
            if ~strcmp(old, vn), S.models = rmfield(S.models, old); end
        else
            error('S.models.%s missing and no case-insensitive match found', vn);
        end
    end
end

% Harmonize classes to reduce memory spikes
likeH = 'double';
if opt.cast_like_H, likeH = class(S.H); end
castlike = @(A) cast(A, likeH);

% ----------- Load & align advection -----------
D = [];
if ~isempty(opt.advection_file) && isfile(opt.advection_file)
    logf('Loading advection map: %s (%s)', opt.advection_file, opt.advection_var);
    A = load(opt.advection_file);
    if isfield(A, opt.advection_var)
        D0 = A.(opt.advection_var);
        if strcmpi(opt.advection_units,'Wm2'), D0 = 1000 * D0; end % ⇒ mW/m^2

        % Try to align to S grid
        if all(isfield(A, {'Xgrid','Ygrid'})) ...
           && isequal(size(A.Xgrid), size(D0)) && isequal(size(A.Ygrid), size(D0))
            D = regrid_to(S.Xgrid, S.Ygrid, A.Xgrid, A.Ygrid, D0);
        elseif isequal(size(D0), size(S.H))
            D = D0;
        elseif isequal(size(D0), [size(S.H,2) size(S.H,1)])
            D = D0.';
        else
            warning('Advection size %dx%d vs model %dx%d. Skipping advection.', ...
                size(D0,1), size(D0,2), size(S.H,1), size(S.H,2));
            D = [];
        end
        if ~isempty(D)
            D = castlike(D);
            if ~isfinite(D(1)), D(~isfinite(D)) = 0; end % safe fill if needed
        end
    else
        warning('Variable "%s" not found in %s. Skipping advection.', opt.advection_var, opt.advection_file);
    end
else
    logf('No advection file provided — computing non-advected fields only.');
end
S.Delta_q_adv = D;

% ----------- Compute Gmin  -----------
assert(exist('processGHF','file')==2, 'processGHF.m not on path.');
logf('Computing Gmin (no-advection)...');
[~, Gmin_noadv, ~] = processGHF(castlike(S.models.Hazzard), castlike(S.Ts), castlike(S.Mb), castlike(S.H));
S.Gmin = castlike(Gmin_noadv);
S.Gmin_adv = [];
if ~isempty(D)
    try
        S.Gmin_adv = castlike(Gmin_noadv) - D;   % q* = q - Δq_adv (mW/m^2)
    catch 
        fprintf('Gmin adv apply failed');
        S.Gmin_adv = [];
    end
end

% ----------- Per-model GΔ -----------
nM = numel(validNames);
if ~isfield(S,'Gdiff'), S.Gdiff = struct(); end
if ~isempty(D) && opt.keep_gdiff_noadv && ~isfield(S,'Gdiff_noadv'), S.Gdiff_noadv = struct(); end

logf('Building GΔ fields for %d models%s ...', nM, tern(~isempty(D), ' (with advection)', ''));
for i = 1:nM
    field_i = validNames{i};
    Mi = castlike(S.models.(field_i));
    [~, ~, Gdif_na] = processGHF(Mi, castlike(S.Ts), castlike(S.Mb), castlike(S.H)); % (M - Gmin_noadv) by your function contract
    Gdif_na = castlike(Gdif_na);

    if ~isempty(D)
        if opt.keep_gdiff_noadv, S.Gdiff_noadv.(field_i) = Gdif_na; end
        % Advected difference = (M - Gmin_noadv) + Δq  ==  M - (Gmin_noadv - Δq)
        S.Gdiff.(field_i) = castlike(Gdif_na + D);
    else
        S.Gdiff.(field_i) = Gdif_na;
    end

    if opt.verbose, fprintf('[add_g_fields] %02d/%02d  %-24s ok\n', i, nM, field_i); end
end

% ----------- Save -----------
outfile = char(opt.outfile); if isempty(outfile), outfile = datafile; end
if opt.save
    logf('Saving to %s (incremental)...', outfile);
    try
        mf = matfile(outfile,'Writable',true);
        mf.S = S;
        logf('Saved -> %s', outfile);
    catch
        fprintf('matfile save failed');
        save(outfile, 'S', '-v7.3');
        logf('Saved (fallback) -> %s', outfile);
    end
end
end

% ---------- helpers ----------
function Dout = regrid_to(Xt, Yt, Xs, Ys, Din)
    % regridding with monotonic axes & class preservation
    D = double(Din);
    Fx = double(Xs(1,:)); Fy = double(Ys(:,1));
    if any(diff(Fx)<=0), Fx = fliplr(Fx); D = fliplr(D); end
    if any(diff(Fy)<=0), Fy = flipud(Fy); D = flipud(D); end
    F  = griddedInterpolant({Fy, Fx}, D, 'linear', 'nearest'); % {row(Y), col(X)}
    Dout = F(double(Yt), double(Xt));
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
