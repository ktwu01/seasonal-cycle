function makeArgoMask(cesmFname, argoFname, maskSaveName)
 
% matlab function to identify contiguous regions in the ocean that are not measured by Argo

% read in heat content for each to identify mask
argoHeat = ncread(argoFname, 'heatContentMap');
argoLat = ncread(argoFname, 'LATITUDE');
argoLon = ncread(argoFname, 'LONGITUDE');
   
cesmHeat = ncread(cesmFname, 'heatContentMap');
cesmLat = ncread(cesmFname, 'lat');
cesmLon = ncread(cesmFname, 'lon');
 
% check that lat and lon are the same
if any(argoLat ~= cesmLat)
	error('grids not aligned in latitude')
elseif any(argoLon ~= cesmLon)
	error('grids not aligned in longitude')
end

mask = isnan(argoHeat(:, :, 1)) & ~isnan(cesmHeat(:, :, 1));

[B L] = bwboundaries(mask > 0); % gives contiguous regions
L(L == 0) = NaN;
% remove values that are 'holes'
L(mask == 0) = NaN;

% rename entries of L to be in descending size order
Lnew = L;
values = unique(Lnew(:));
values(isnan(values)) = [];
instances = histc(Lnew(:),values);
[v idx] = sort(instances, 'descend');
for ct = 1:length(idx)
	Lnew(L == idx(ct)) = ct;
end

% KAM, an issue may emerge because the first cluster spans many latitudes
% may want to fix later

% save to ascii filed
dlmwrite(maskSaveName, mask)

quit




