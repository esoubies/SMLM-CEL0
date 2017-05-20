function normai= ComputeNorm_ai(PSF,coefEch)
% function normai= ComputeNorm_ai(PSF,coefEch)
%
% Computes the norm of the columns of the operator 
% MechT*H*x*Mech where H is the convolution operator 
% defined by the PSF
%
% Copyright (C) 2017 S. Gazagnes sgazagnes@gmail.com, E. Soubies esoubies@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

SRSize = size(PSF,1);
Mech = zeros(SRSize/coefEch, SRSize);
for i = 0 : SRSize/coefEch-1
    Mech(i+1, 1+coefEch*i : coefEch+coefEch*i) = 1;
end
MechT = Mech';

for i = 1:coefEch
    for j = 1:coefEch
        matr = Mech * circshift(PSF, [j-1,i-1]) * MechT;
        normai(i,j) = sqrt(sum(matr(:).^2));
    end
end
normai = repmat(normai, size(Mech,1));
end

