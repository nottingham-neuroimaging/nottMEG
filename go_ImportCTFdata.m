function [resource, data] = go_ImportCTFdata(varargin)
%
% CTF Dataset Import Wrapper - George O'Neill 2015
%
% [resource,data] = go_ImportCTFData(ds,chanType)
%
% Inputs:
%
%   - ds: A string with the path and name of the CTF .ds file.
%   - chanType: {'M'} | 'E' | 'U' | 'HLC'
%               Selecting 'M' (or not specifying) imports MEG data, 'E'
%               will import EEG Channels and 'U' will import ADC/DAC and
%               parallel port triggers. 'HLC' will extrat head localisation
%               information.
%
% Outputs:
%
%   - resource: A header which contains metadata.
%   - data: The CTF data, output in a 2D array.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Importing CTF Data: ')

switch nargin
    case 1
        ds = varargin{1};
        chanType = 'M';
    otherwise
        ds = varargin{1};
        chanType = upper(varargin{2}); 
end

if ~exist(ds)
    error('ERROR: .ds file not found.')
end

hdr = readCTFds(ds);
channels = find(sum(hdr.res4.chanNames(:,1:length(chanType)) == repmat(chanType,length(hdr.res4.chanNames),1),2) == length(chanType));

data = getCTFdata(hdr,[],channels,'fT');

if nargin < 3
    data = permute(data,[1 3 2]);
    data = reshape(data,size(data,1)*size(data,2),size(data,3))';
end

resource = hdr.res4;
fprintf('     COMPLETE\n')
