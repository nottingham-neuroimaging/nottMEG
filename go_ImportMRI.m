function [mri_info, mri_image] = go_ImportMRI(mri)
%
% CTF MRI Import Wrapper - George O'Neill 2015
%
% [mri_info,mri_image] = go_readCTFmri(mri)
%
% Inputs:
%  
%   - mri: A string with the path and name of the CTF .mri file.
%
% Outputs:
%
%   - mri_info: A small header which contains the MRI resolution and
%              transformation matricies required for source reconstruction.
%   - mri_image: The anatomical image, presented in FSL orientation.
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
fprintf('Loading CTF MRI: ')

if ~exist(mri)
    error('ERROR: .mri files does not exist')
end

[MRItag,MRIdata]=readCTFMRI(mri);


for ind = 1:length(MRItag);
    if strmatch(MRItag(ind).name,'_CTFMRI_MMPERPIXEL');
        tmp = MRItag(ind).data;
        tmp(strfind(tmp,'\')) = ';';
        mri_info.VoxSize = str2num(tmp)';
    elseif strmatch(MRItag(ind).name,'_CTFMRI_TRANSFORMMATRIX');
        tmp = MRItag(ind).data;
        tmp(strfind(tmp,'\')) = ';';
        mri_info.T = reshape(str2num(tmp),4,4);
    end      
end


anatomy = permute(MRIdata,[3 1 2]);
anatomy = flipdim(anatomy,1);
anatomy = flipdim(anatomy,2);
anatomy = flipdim(anatomy,3);

mri_image = anatomy;

fprintf('        COMPLETE\n')

