function [weights, leads] = go_BeamformerFixedOrientation(hdm,inv_Cr,r,or)
%
% Beamformer of a dipole of known orientation - George O'Neill 2016
%
% Beamformer assumes a localsphere forward model is being used, generated
% from go_ft_Localspheres. Lead field calculation is based on Sarvas, PMB
% 1987, and written in MEX-complied C for speeeeeeed.
%
% [weights, leads] = go_BeamformerFixedOrientation(hdm,inv_Cr,r,or)
%
% Inputs:
%
%   - hdm: Head model structure which is generated in go_ft_Localspheres.
%   - inv_Cr: Inverted regularised covariance matrix used in beamfomer
%             weights generation.
%   - r: 1x3 vector containing the location of interest (must be in cm and 
%        CTF space)
%   - or: 1x3 vector with the dipole orientation, vector is automatically
%         scaled to be a unit vector in not already.
%
% Outputs:
%
%   - weights: Vector containing the beamformer weights for reconstruction
%   - leads: Vector containing associated lead fields
%
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

% Check if r is a three vector
if length(r) ~= 3
    error('Dipole location must be a 1x3 vector input')
else
    dippos = r;
end

if length(or) ~= 3
    error('Dipole orientation must be a 1x3 vector input')
end

%Unpack headmodel structure for clarity.
sens = hdm.sens;

% Calculate the lead fields for three orthogonally oriented dipoles.

ncoils = length(sens.coilpos);

if size(hdm.r, 1)~=ncoils
    error('number of spheres is not equal to the number of coils')
end

if size(hdm.o, 1)~=ncoils
    error('number of spheres is not equal to the number of coils');
end

% normalise the orientation if not a unit vector
or = or./norm(or);

lf = zeros(ncoils,1);
% Use the multiple spheres forward solution
for coil=1:ncoils
    for dip=1
        % shift dipole and magnetometer coil to origin of sphere
        tmppos  = dippos(dip, :) - hdm.o(coil, :);
        coilpos = sens.coilpos(coil, :) - hdm.o(coil, :);
        lf(coil) = go_leadfield(tmppos, or, coilpos, sens.coilori(coil, :));
    end
end


if isfield(sens, 'tra')
    % this appears to be the modern complex gradiometer definition
    % construct the channels from a linear combination of all magnetometers
    lf = sens.tra * lf;
    
end

leads = lf(sens.chanidx);
weights = (leads'*inv_Cr/(leads'*inv_Cr*leads))';
