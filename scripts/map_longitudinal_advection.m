function map_longitudinal_advection()
% Builds Δq_adv map (W/m^2) with robust masking + validation
%
% Inputs expected:
%   datasets/Ts_interp.mat         -> Ts_interp (K)
%   datasets/Mb_interp.mat         -> Mb_interp (m ice eq / yr)
%   datasets/coldex_icethk.mat     -> Xgrid (m), Ygrid (m), H (m)
%   datasets/spb_merged.shp        -> South Pole Basin polygons (EPSG:3031)
%   datasets/mouginot_icevel.mat   -> VX, VY on (x,y) in EPSG:3031 meters, m/yr units
%
% Output:
%   datasets/longitudinal_advection_maps.mat
%   datasets/longitudinal_advection_SPB.mat

%% -------------------- Settings --------------------
clearvars -except varargin; clc;
in_dir      = 'datasets';
f_no_sliding  = 0.5;               % 0.5 (shear-dominated) or 1.0 (plug/sliding)
outdir     = 'datasets';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% -------------------- Load --------------------
fprintf('Loading core inputs...\n');
load(fullfile(in_dir,'coldex_icethk.mat'),'Xgrid','Ygrid','H');  % meters
load(fullfile(in_dir,'Ts_interp.mat'),'Ts_interp');              % K
load(fullfile(in_dir,'Mb_interp.mat'),'Mb_interp');              % m ice eq / yr
Ts   = Ts_interp;
adot_yr = Mb_interp;                                             % m/yr 
dx = mean(diff(Xgrid(1,:)));
dy = mean(diff(Ygrid(:,1)));
sec_per_yr = 365.25*24*3600;

V = load(fullfile(in_dir,'mouginot_icevel.mat'));  

% Ensure grids match
assert(isequal(size(V.Xgrid),size(Xgrid)) && isequal(size(V.Ygrid),size(Ygrid)), ...
    'Velocity grid size does not match Xgrid/Ygrid.');

% Pull components if present
have_components = isfield(V,'vx') && isfield(V,'vy') && ~isempty(V.vx) && ~isempty(V.vy);

if have_components
    vx_yr = V.vx;              % m/yr on your grid
    vy_yr = V.vy;              % m/yr
    U_yr  = hypot(vx_yr, vy_yr);
    % Unit flow vectors from components
    Umag = max(U_yr, 1e-9);
    ex = vx_yr ./ Umag;
    ey = vy_yr ./ Umag;
else
    bearing_from_north_deg = 25;            
    th = deg2rad(bearing_from_north_deg);
    ex = sin(th) * ones(size(Xgrid));
    ey = cos(th) * ones(size(Xgrid));

    assert(isfield(V,'icevel'), 'Need icevel (speed) if vx/vy are absent.');
    U_yr = V.icevel;           % m/yr
end

% Convert speed to m/s and pick depth-mean factor
U_s   = U_yr / sec_per_yr;     % m/s
u_bar = f_no_sliding .* U_s;   % m/s

%% -------------------- Directional derivatives --------------------
[dTs_dx, dTs_dy] = gradient(Ts, dx, dy);
[dH_dx,  dH_dy ] = gradient(H,  dx, dy);
[da_dx,  da_dy ] = gradient(adot_yr, dx, dy);

dTs_along = dTs_dx .* ex + dTs_dy .* ey;      % K/m
dH_along  = dH_dx  .* ex + dH_dy  .* ey;      % unitless
da_along  = da_dx  .* ex + da_dy  .* ey;      % (m/yr)/m

adot = adot_yr / sec_per_yr;                  % m/s
da_along = da_along / sec_per_yr;             % (m/s)/m

%% -------------------- Material props & basal gradient G0 --------------------
rho   = 917;                                  % kg/m^3
g     = 9.81;                                 % m/s^2
gamma = -7.42e-8;                             % K/Pa (melting point depression)
c     = 152.5 + 7.122*Ts;                     % J/kg/K  (evaluate at Ts)

T_melt_base = 273.15 + gamma * rho * g .* H;  % K  (pressure-melting estimate)
G0 = (T_melt_base - Ts) ./ max(H,1);          % K/m (first-pass)
G0(~isfinite(G0)) = 0; G0 = max(G0,0);

%% -------------------- Longitudinal advection equivalent basal flux --------------------
term_geom = ( (1./max(H,1)).*dH_along ) - ( (1./max(adot,1e-12)).*da_along );
Delta_q_adv = rho .* c .* u_bar .* ( H .* dTs_along + 0.333 .* G0 .* H.^2 .* term_geom ); % W/m^2

% Mask invalid 
bad = ~isfinite(Delta_q_adv) | ~isfinite(dTs_along) | ~isfinite(u_bar) | (H<=0) | ~isfinite(adot);
Delta_q_adv(bad) = NaN;
Delta_q_adv_mW = 1e3 * Delta_q_adv;

%% -------------------- Maps --------------------
figure('Color','w','Position',[80 80 1100 480]);

imagesc(Xgrid(1,:)/1000, Ygrid(:,1)/1000, Delta_q_adv_mW);
axis image; set(gca,'YDir','normal');
c = colorbar;
ylabel(c, '\Delta q_{adv} (mW m^{-2})');
xlabel('X (km)'); ylabel('Y (km)');
title('\Delta q_{adv} (mW m^{-2})');

%% -------------------- Stats + central-95% histograms --------------------
vals_all = Delta_q_adv(:); vals_all = vals_all(isfinite(vals_all));
if isempty(vals_all)
    error('No finite Δq_{adv} after masking—check inputs/masks.');
end

p_full = prctile(vals_all,[1 5 25 50 75 95 99]);
mu = mean(vals_all); sg = std(vals_all);

fprintf('\nΔq_adv summary (W/m^2):\n');
fprintf('  min/1%%/5%%/25%%/50%%/75%%/95%%/99%%/max =\n');
fprintf('  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n', ...
        min(vals_all), p_full(1), p_full(2), p_full(3), p_full(4), ...
        p_full(5), p_full(6), p_full(7), max(vals_all));
fprintf('  mean = %.4g, std = %.4g\n', mu, sg);

p95 = prctile(vals_all,[2.5 97.5]);
vals_clip    = vals_all(vals_all>=p95(1) & vals_all<=p95(2));
vals_mW_clip = vals_clip * 1e3;

histogram(vals_mW_clip, 80, 'Normalization','pdf');
xlabel('\Delta q_{adv} (mW m^{-2})'); ylabel('PDF');
title('Central 95% \Delta q_{adv} (mW m^{-2})'); grid on;

%% -------------------- Save outputs --------------------
save(fullfile(outdir,'longitudinal_advection_maps.mat'), ...
     'Delta_q_adv','dTs_along','Xgrid','Ygrid','H','u_bar','G0','ex','ey');

write_xyz(fullfile(outdir,'Delta_q_adv.xyz'), Xgrid, Ygrid, Delta_q_adv);
write_xyz(fullfile(outdir,'dTs_along.xyz'),   Xgrid, Ygrid, dTs_along);

fprintf('\nDone. Outputs in %s\n', outdir);

end

%  helper
function write_xyz(fname, X, Y, Z)
    fid = fopen(fname,'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    Xv = X(:); Yv = Y(:); Zv = Z(:);
    good = isfinite(Zv);
    fprintf(fid,'%.3f %.3f %.6g\n',[Xv(good) Yv(good) Zv(good)]');
    fclose(fid);
end
