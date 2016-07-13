function hdm = go_ft_LocalSpheres(ds)

% Huang Head Model reader - George O'Neill 2016
% Uses the FieldTrip distribution from 16 Jan 2016.
% hdm = go_SingleShellModel(ds);
%
% Inputs:
%
%     -ds: Path to the CTF dataset
%
% Outputs:
%
%   - hdm: Structure which contains headmodel information
%              transformation matricies required for source reconstruction
%              including sensor information.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slsh = strfind(ds,'/');
base = ds(1:slsh(end));

% Check to see if the head model exists and load, otherwise let
% fieldtrip generate model and forward parameters.

if exist([ds '/localspheres.mat']) == 0
    %     warning('OFF')
    addpath '/net/geraint/data_local/george/fieldtrip/fieldtrip-master/';
    ft_defaults
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load gradiometer and shape information from data
    
    hdr = ft_read_header(ds);
    if ~exist([ds '/default.hdm'])
        error('Local spheres head model not generated, please use LocalSpheres in command line')
    else
        headmodel = ft_read_vol([ds '/default.hdm']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % headmodel
    sens = ft_read_sens(ds);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If the volume conduction model consists of multiple spheres then we
    % have to match the channels in the gradiometer array and the volume
    % conduction model.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the initial localspheres volume conductor has a local sphere per
    % channel, whereas it should have a local sphere for each coil
    if size(headmodel.r,1)==size(sens.coilpos,1) && ~isfield(headmodel, 'label')
        % it appears that each coil already has a sphere, which suggests
        % that the volume conductor already has been prepared to match the
        % sensor array
        return
    elseif size(headmodel.r,1)==size(sens.coilpos,1) && isfield(headmodel, 'label')
        if ~isequal(headmodel.label(:), sens.label(:))
            % if only the order is different, it would be possible to reorder them
            error('the coils in the volume conduction model do not correspond to the sensor array');
        else
            % the coil-specific spheres in the volume conductor should not have a label
            % because the label is already specified for the coils in the
            % sensor array
            headmodel = rmfield(headmodel, 'label');
        end
        return
    end
    
    % the CTF way of representing the headmodel is one-sphere-per-channel
    % whereas the FieldTrip way of doing the forward computation is one-sphere-per-coil
    Nchans   = size(sens.tra,1);
    Ncoils   = size(sens.tra,2);
    Nspheres = size(headmodel.label);
    
    if isfield(headmodel, 'orig')
        % these are present in a CTF *.hdm file
        singlesphere.o(1,1) = headmodel.orig.MEG_Sphere.ORIGIN_X;
        singlesphere.o(1,2) = headmodel.orig.MEG_Sphere.ORIGIN_Y;
        singlesphere.o(1,3) = headmodel.orig.MEG_Sphere.ORIGIN_Z;
        singlesphere.r      = headmodel.orig.MEG_Sphere.RADIUS;
        % ensure consistent units
        singlesphere = ft_convert_units(singlesphere, headmodel.unit);
        % determine the channels that do not have a corresponding sphere
        % and use the globally fitted single sphere for those
        missing = setdiff(sens.label, headmodel.label);
        if ~isempty(missing)
            warning('using the global fitted single sphere for %d channels that do not have a local sphere', length(missing));
        end
        for i=1:length(missing)
            headmodel.label(end+1) = missing(i);
            headmodel.r(end+1,:)   = singlesphere.r;
            headmodel.o(end+1,:)   = singlesphere.o;
        end
    end
    
    % make a new structure that only holds the local spheres, one per coil
    localspheres = [];
    localspheres.type = 'localspheres';
    localspheres.unit = headmodel.unit;
    
    % for each coil in the MEG helmet, determine the corresponding channel and from that the corresponding local sphere
    for i=1:Ncoils
        coilindex = find(sens.tra(:,i)~=0); % to which channel does this coil belong
        if length(coilindex)>1
            % this indicates that there are multiple channels to which this coil contributes,
            % which happens if the sensor array represents a synthetic higher-order gradient.
            [dum, coilindex] = max(abs(sens.tra(:,i)));
        end
        
        coillabel = sens.label{coilindex};               % what is the label of this channel
        chanindex = find(strcmp(coillabel, headmodel.label));  % what is the index of this channel in the list of local spheres
        localspheres.r(i,:) = headmodel.r(chanindex);
        localspheres.o(i,:) = headmodel.o(chanindex,:);
    end
    hdm = localspheres;
    sens.chanidx = strmatch('meggrad',sens.chantype);
    hdm.sens = sens;
    
    % save headmodel in CTF folder.
    try
    save([ds '/localspheres.mat'],'hdm')
    catch
        sprintf('WARNING: extracted headmodel cannot be saved to CTF folder\n')
    end
    warning('ON')
    
else
    load([ds '/localspheres.mat'])
end
