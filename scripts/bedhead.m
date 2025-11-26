function varargout = bedhead(varargin)
% bedhead returns static pressure head at an ice sheet base in units of meters 
% freshwater equivalent.  This function calculates only a simple first-order 
% approximation and does not consider effects such as velocity or strain or 
% really anything fancy. 
% 
%% Requirements 
% This function requires no special toolboxes and can be used for any glacier
% or ice sheet DEM, but if you do not supply your own DEM, the bedhead function
% will automatically use the Bedmap2 DEM.  If you do not supply a DEM of your 
% own you'll need the Bedmap2 Toolbox for Matlab and Antarctic Mapping Tools,
% both of which are available on the Mathworks File Exchange site.   
% 
%% Citation
% There is only one simple mathematical expression in this function, which I adapted
% from Section 3.3 of 
% 
% Fricker, H.A. and T. Scambos. "Connected subglacial lake activity on lower Mercer
% and Whillans ice streams, West Antarctica, 2003?2008." Journal of Glaciology 55.190 
% (2009): 303-315,
% 
% Fricker and Scambos and just about everyone else who has done work in this area cite 
% Shreve, who seems to be the originator of the model: 
% 
% Shreve, R. L. "Movement of water in glaciers." Journal of Glaciology 11 (1972): 205-214.
% 
%% Syntax 
% 
%  psi = bedhead
%  psi = bedhead(lat,lon)
%  psi = bedhead(x,y) 
%  psi = bedhead(...,extrakm) 
%  psi = bedhead(...,'surface',SurfaceDEM)
%  psi = bedhead(...,'bed',BedDEM)
%  psi = bedhead(...,'thickness',ThicknessDEM)
%  psi = bedhead(...,'rhoi',IceDensity)
%  psi = bedhead(...,'rhow',WaterDensity)
%  [X,Y,psi] = bedhead(...)
% 
%% Description 
% 
% psi = bedhead returns 6667x6667 1 km gridded static pressure head for Antarctica from Bedmap2. 
%
% psi = bedhead(lat,lon) returns 1 km gridded static pressure head from Bedmap2 bounded by a polar
% stereographic quadrangle just large enough to encompass all the points in lat,lon. This syntax 
% follows the bedmap2_data function. 
% 
% psi = bedhead(x,y)  returns 1 km gridded static pressure head from Bedmap2 bounded by a polar
% stereographic quadrangle just large enough to encompass all the polar stereographic (ps71 meters) 
% coordinates x,y. This syntax follows the bedmap2_data function. 
% 
% psi = bedhead(...,extrakm) adds a specified extra number of kilometers around x,y or lat,lon if the
% syntax above is used.  This adds a buffer around your data, or if you would like a 500-km wide grid 
% centered on a single point (75°S,100°E) let extrakm equal 250, i.e., bedhead(-75,100,250), 
% 
% psi = bedhead(...,'surface',SurfaceDEM) specifies a gridded surface elevation DEM. If you specify a 
% surface DEM you must also specify a BedDEM or a ThicknessDEM. 
% 
% psi = bedhead(...,'bed',BedDEM) specifies a gridded bed elevation DEM. If you specify a bed elevation
% DEM you must also specify a SurfaceDEM or a ThicknessDEM. 
% 
% psi = bedhead(...,'thickness',ThicknessDEM) specifies a gridded ice thickness DEM. If you specify a 
% thickness DEM you must also specify a BedDEM or a SurfaceDEM. 
% 
% psi = bedhead(...,'rhoi',IceDensity) specifies ice density. If you specify your own DEM you may also
% define IceDensity as a spatially-variable grid, say if you want to fully model variations in column-
% averaged density due to firn air content. Default value is 917 kg/m^3. 
% 
% psi = bedhead(...,'rhow',WaterDensity) specifies density of fresh water. Default value is 1000 kg/m^3. 
% 
% [X,Y,psi] = bedhead(...) returns gridded polar stereographic meters X, Y corresponding to the psi grid.  
% 
%% Example 1: Hydrostatic potential under Kamb Ice Stream: 
% Get a 600-km wide map of hydrostatic potential under Kamb Ice Stream: 
% 
%   [kamblat,kamblon] = scarloc('kamb ice stream');
%   [X,Y,psi] = bedhead(kamblat,kamblon,300); 
% 
%   pcolor(X,Y,psi)
%   shading flat
%   axis equal
%   scarlabel('Kamb Ice Stream','fontangle','italic') 
%   cb = colorbar; 
%   ylabel(cb,'hydraulic head (m)') 
% 
%% Example 2: Define your own DEM:
% Suppose you have some arbitrary 1 km grid, 400 km wide, centered on Byrd Camp.   
% 
%   [lat,lon] = psgrid('byrd camp',400,1); 
% 
% To make this example easy I'm using Bedmap2 surface and bed elevations as a starting point: 
% 
%   sfz = bedmap2_interp(lat,lon,'surface'); 
%   bed = bedmap2_interp(lat,lon,'bed'); 
% 
% But Bedmap2 has its uncertainties. Let's say it's 10 m for the surface and the bed 
% uncertainty is given by: 
% 
%   bedunc = bedmap2_interp(lat,lon,'beduncertainty'); 
% 
% Create surface and surface and bed DEMs with random values added to simulate 
% uncertainties (crude approach, I know): 
% 
%   sfz2 = sfz + 10*randn(size(sfz)); 
%   bed2 = bed + bedunc.*randn(size(bed)); 
% 
% Pressure head is then given by: 
% 
%   psi2 = bedhead('surface',sfz2,'bed',bed2); 
% 
%% Author Info
% This function and supporting documentation were written by Chad A. Greene of the University
% of Texas at Austin's Institute for Geophysics (UTIG), January 2016.  I did not come up with 
% the model for subglacial hydraulic head. If you use this function please cite the Shreve 1972 
% paper listed above. 
% http://www.chadagreene.com

%% Set defaults: 

rhoi = 917;   % ice density (kg/m^3)
rhow = 1000;  % water density (kg/m^3)
usebedmap2 = true; 
extrakm = 0; 

%% Parse inputs: 

if nargin>0 
    
    % Subset the dataset by region?  Assume regional subsetting if first input is numeric:   
    if isnumeric(varargin{1})
        subsetregion = true; 
        assert(isnumeric(varargin{2})==1,'Input error: If hydropotential''s first input is numeric I assume the first two inputs are coordinates. But you have entered non-numeric inputs for the second coordinate.') 
        lati_or_xi = varargin{1}; 
        loni_or_yi = varargin{2}; 
        if nargin>2
            if isnumeric(varargin{3})
                extrakm = varargin{3}; 
            end
        end
    end
    
    % Ice density: 
    tmp = strcmpi(varargin,'rhoi'); 
    if any(tmp)
        rhoi = varargin{find(tmp)+1}; 
    end
    
    % Water density: 
    tmp = strcmpi(varargin,'rhow'); 
    if any(tmp)
        rhow = varargin{find(tmp)+1}; 
    end
    
    % Surface elevation: 
    tmp = strncmpi(varargin,'surface',3); 
    if any(tmp)
        sfz = varargin{find(tmp)+1}; 
        usebedmap2 = false; 
    end
    
    % Bed elevation: 
    tmp = strcmpi(varargin,'bed'); 
    if any(tmp)
        bed = varargin{find(tmp)+1}; 
        usebedmap2 = false; 
    end
    
    % ice thickness: 
    tmp = strncmpi(varargin,'thickness',3); 
    if any(tmp)
        thickness = varargin{find(tmp)+1}; 
        usebedmap2 = false; 
    end
        
end

%% Build or import DEMs: 


if usebedmap2 
    if subsetregion
        [X,Y,bed] = bedmap2_data('bed',lati_or_xi,loni_or_yi,extrakm,'xy'); 
        thickness = bedmap2_data('thickness',lati_or_xi,loni_or_yi,extrakm,'xy'); 
        mask = bedmap2_data('icemask',lati_or_xi,loni_or_yi,extrakm,'xy'); 
        
    else
        [X,Y,bed] = bedmap2_data('bed','xy'); 
        thickness = bedmap2_data('thickness','xy'); 
        mask = bedmap2_data('icemask','xy'); 
    end
else
    
    assert(exist('bed','var')+exist('sfz','var')+exist('thickness','var')>1,'Input error: DEMs are not fully defined or are of mixed dimensions.') 

    if ~exist('thickness','var')
        thickness = sfz-bed; 
    end

    if ~exist('bed','var') 
        bed = sfz-thickness; 
    end
end

%% Perform mathematics: 

head = thickness.*rhoi./rhow + bed; 

% Ensure ice shelves and open ocean are NaN:   
if usebedmap2
    head(mask~=0) = NaN; 
end

%% Define function outputs: 

switch nargout
    case 1
        varargout{1} = head; 
    case 3
        varargout{1} = X; 
        varargout{2} = Y; 
        varargout{3} = head; 
    otherwise
        error('unrecognized number of outputs') 
end


end
