function make_sink_mask()
% Build sink mask from hydraulic head on a 2 m grid, then
% regrid to 100 m
%
% Outputs:
%   datasets/sink_mask_new.mat  -> sink_mask, Xgrid, Ygrid
%   datasets/sink_mask_comp.mat -> sink_mask_comp, comp_id, Xgrid, Ygrid

%% -------------------- CONFIG --------------------
conn        = 8;                      % 4 or 8 connectivity (for 100 m comps)
min_size    = 5;                      % min pixels per component
out_sinks   = fullfile('datasets_for_gmin','sink_mask_new.mat');
out_comp    = fullfile('datasets_for_gmin','sink_mask_comp.mat');


%% -------------------- Load --------------------
disp('Loading 2 m input data...');
load('datasets/coldex_bedelv_2m.mat',  'B_2m', 'X2m', 'Y2m');
load('datasets/REMA_srfelv_2m.mat',    'S_rem_2m');
load('datasets/coldex_srfelv_2m.mat',  'S_coldex_2m');

assert(isequal(size(B_2m), size(S_rem_2m), size(S_coldex_2m), size(X2m), size(Y2m)), ...
    'Size mismatch among 2 m rasters');

% Basic diagnostics
xv2 = X2m(1,:);  
yv2 = Y2m(:,1);
dx2 = median(diff(xv2)); 
dy2 = median(diff(yv2));
fprintf('[2m grid] size = %dx%d | X %s (dx=%.3f m) | Y %s (dy=%.3f m)\n', ...
    size(B_2m,1), size(B_2m,2), tern(diff(xv2)>0,'asc','desc'), dx2, tern(diff(yv2)>0,'asc','desc'), dy2);

%% -------------------- Ice thickness --------------------
proj = projcrs(3031); 
[lat2, ~] = projinv(proj, X2m, Y2m);

mask_coldex = (lat2 <= -88);    % use COLDEX surface near the pole, REMA elsewhere
combined_srfelv_2m = S_rem_2m;
combined_srfelv_2m(mask_coldex) = S_coldex_2m(mask_coldex);

% ---- Area where COLDEX replaces REMA (south of 88°S) ----
valid_coldex = mask_coldex & isfinite(S_coldex_2m);

pixel_area_m2   = abs(dx2 * dy2);
n_coldex_pixels = nnz(valid_coldex);
area_coldex_km2 = n_coldex_pixels * pixel_area_m2 / 1e6;

n_finite = nnz(isfinite(combined_srfelv_2m));
frac_coldex_pct = 100 * n_coldex_pixels / n_finite;

fprintf(['[COLDEX replace] %d pixels -> %.2f km^2 ' ...
         '(%.2f%% of finite 2 m surface grid)\n'], ...
        n_coldex_pixels, area_coldex_km2, frac_coldex_pct);

% Ice thickness
thk_2m = max(combined_srfelv_2m - B_2m, 0);

%% -------------------- Normalize to north-up for GRIDobj --------
flipY2 = yv2(1) < yv2(end);  
if flipY2
    disp('[2m orient] Detected south-up input (Y ascending). Flipping to north-up for processing.');
    Bz2     = flipud(B_2m);
    Thkz2   = flipud(thk_2m);
    yv2_north = flipud(yv2);    
else
    Bz2     = B_2m;
    Thkz2   = thk_2m;
    yv2_north = yv2;            
end
xv2_east = xv2;                 

%% -------------------- GRIDobjs for TopoToolbox -----------------
disp('Creating 2 m GRID objects (north-up) ...');
bed_dem_2m = GRIDobj(xv2_east, yv2_north, Bz2);
thk_dem_2m = GRIDobj(xv2_east, yv2_north, Thkz2);

%% -------------------- Hydraulic head ---------------------
try
    disp('Calculating hydraulic head with bedhead() on 2 m grid ...');
    head_dem_2m = bedhead('bed', bed_dem_2m, 'thickness', thk_dem_2m);
catch
    warning('bedhead() not found; using bed elevation as head proxy (2 m).');
    head_dem_2m = bed_dem_2m;
end

%% -------------------- Fill sinks -------------------------
disp('Computing DEM fill (2 m, north-up) ...');
DEMf_2m = fillsinks(head_dem_2m);

sinks_depth_north_2m = DEMf_2m.Z - head_dem_2m.Z;
sink_mask_north_2m   = (sinks_depth_north_2m > 0) & isfinite(sinks_depth_north_2m);

if flipY2
    sink_mask_2m = flipud(sink_mask_north_2m);
else
    sink_mask_2m = sink_mask_north_2m;
end

%% -------------------- Regrid sinks to 100 m analysis grid ------------
% The bootstrap / GHF evaluation uses S.Xgrid,S.Ygrid (100 m)
disp('Regridding 2 m sink mask to 100 m analysis grid ...');

L100 = load('datasets_for_gmin/coldex_bedelv.mat', 'B', 'Xgrid', 'Ygrid');
Xgrid = L100.Xgrid;
Ygrid = L100.Ygrid;

xv2a = xv2;
yv2a = yv2;
sink2a = sink_mask_2m;
if numel(xv2a)>1 && xv2a(2) < xv2a(1)
    xv2a  = fliplr(xv2a);
    sink2a = fliplr(sink2a);
end
if numel(yv2a)>1 && yv2a(2) < yv2a(1)
    yv2a  = flipud(yv2a);
    sink2a = flipud(sink2a);
end

% Interpolate
F_sink = griddedInterpolant({yv2a, xv2a}, double(sink2a), 'linear', 'none');
sink_interp = F_sink(Ygrid, Xgrid);
sink_interp(~isfinite(sink_interp)) = 0;

sink_mask = sink_interp > 0.5;   

fprintf('[100m sinks] grid size = %dx%d | sink px = %d (%.3f%%)\n', ...
    size(sink_mask,1), size(sink_mask,2), nnz(sink_mask), ...
    100*nnz(sink_mask)/numel(sink_mask));

%% -------------------- Save --------------------
if ~exist(fileparts(out_sinks), 'dir'); mkdir(fileparts(out_sinks)); end
save(out_sinks, 'sink_mask', 'Xgrid', 'Ygrid', '-v7.3');
fprintf('Saved 100 m sink mask to: %s\n', out_sinks);

%% -------------------- Connected components on 100 m grid -------------
disp('Finding connected components on 100 m grid ...');
CC = bwconncomp(sink_mask, conn);
fprintf('  initial components: %d\n', CC.NumObjects);

keepC = cellfun(@numel, CC.PixelIdxList) >= min_size;
sink_mask_comp = false(size(sink_mask));
if any(keepC)
    sink_mask_comp(vertcat(CC.PixelIdxList{keepC})) = true;
end

CC2 = bwconncomp(sink_mask_comp, conn);
comp_id = uint32(labelmatrix(CC2));
fprintf('  kept components: %d (removed %d)\n', CC2.NumObjects, CC.NumObjects - CC2.NumObjects);

%% -------------------- Save component labels -------------------
if ~exist(fileparts(out_comp), 'dir'); mkdir(fileparts(out_comp)); end
save(out_comp, 'sink_mask_comp', 'comp_id', 'Xgrid', 'Ygrid', '-v7');
fprintf('Saved 100 m component mask/labels -> %s\n', out_comp);

%% -------------------- Quicklook ------------
figure('Color','w','Name','Kept components over 100 m bed');
imagesc(Xgrid(1,:), Ygrid(:,1), L100.B); 
axis image; set(gca,'YDir','normal'); colormap gray; hold on
title(sprintf('Kept sink components (100 m; n=%d)', CC2.NumObjects));
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on
contour(Xgrid, Ygrid, double(sink_mask_comp), [0.5 0.5], 'm', 'LineWidth', 2);

end

% --- helper ---
function s = tern(cond, a, b)
if all(cond), s = a; else, s = b; end
end
