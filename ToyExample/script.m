%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- TOY EXAMPLE SMLM-CEL0
%
% Please read the function SMLMCEL0 instructions to 
% see and understand the options. 
%
% Copyright (C) 2017 S. Gazagnes sgazagnes@gmail.com, 
%                    E. Soubies esoubies@gmail.com
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all , close all;
addpath('../')
% === +++++++ === %
% === Options === %
% === +++++++ === %
options.coefEch = 4;         % Data pixel subdivision coefficient 
options.x0 = 0;              % Initial point for optimization (all 0)
options.lambda = 3.5e-2;     % Regularization parameters to use
options.itmaxIRL=200;        % Max iterations outer loop (IRL1)
options.itmaxFista = 1000;   % Max iterations inner loop (FISTA)
options.verbose=1;           % Verbose mode
options.costFunction = 0;    % Do not compute the cost a each iteration          
psf=load('PSF'); options.PSF=psf.pp/sum(psf.pp(:));          % Set PSF
options.normai=ComputeNorm_ai(options.PSF,options.coefEch);  % Norm of the columns of the operator

% === ++++++++++++++ === %
% === Reconstruction === %
% === ++++++++++++++ === %
Acq = fitsread('NoisyAcquisition.fits'); Acq=Acq/max(Acq(:));   % Read Data
[Recons, infos] = SMLMCEL0 (Acq,options);                       % Run algorithm

% === ++++++++++++++++++++++++++++++ === %
% ===   Comparison with Groud Truth  === %
% === ++++++++++++++++++++++++++++++ === %
gt = fitsread('GroundTruth.fits');
[Y,X]=ind2sub([size(gt,1),size(gt,2)],find(gt>0));
figure;
NbR=length(Recons);
for i=1:NbR
    [Yr,Xr]=ind2sub([size(Recons{i},1),size(Recons{i},2)],find(Recons{i}>0));
    subplot(1,NbR,i);imagesc(Acq);
    axis image; axis off;hold all;colormap gray;
    plot((Xr+0.5)/options.coefEch,(Yr+0.5)/options.coefEch,'xr');
    plot((X+0.5)/options.coefEch,(Y+0.5)/options.coefEch,'og');
    title(['Recons (Red) GT (Green) \lambda=',num2str(options.lambda(i))]);
end
