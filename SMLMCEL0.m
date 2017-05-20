function [SRX, infos] = SMLMCEL0(Iacq, options)
% function [SRX, infos] = SMLMCEL0(Iacq, options)
%
% Single Molecule Localization Microscopy [1]: this algorithm allows to
% reconstruct super-resolved images from one acquisition by minimizing:
% 0.5||Ax - y||^2 + \phi_{celo}(x) + i_{>0}(x)
% where A = M*H with M a downsizing operator and H a convolution,
% \phi_{celo} the CELO penalty and i_{>0} the indicator function over
% positive vectors. The minimization is done with the IRL1 algorithm
%
% [1] Simon Gazagnes, Emmanuel Soubies, Laure Blanc-Féraud,
%     "High Density Molecule Localization for Super-Resolution
%     Microscopy Using CEL0 Based Sparse Approximation"  IEEE International
%     Symposium on Biomedical Imaging, 2017.
%
% Inputs - Iacq   : Acquisition to be processed
%        - options: structure with the following fields:
%           -> PSF         : matrix of the Point Spread Function 
%           -> coefEch     : Ratio between acquisition grid and reconstruction
%                            grid. Example: 64x64 pixels acquisitions grid
%                            to be reconstructed on a 256x256 pixels grid:
%                            coefEch has to be equal to 4.
%           -> normai      : norm of the columns of the linear operator A
%                            if empty, will be computed with the function 
%                            ComputeNorm_ai.m
%           -> lambda      : vector of lambdas (regularization parameter)
%                            to be tested.
%           -> x0          : starting point of the algorithm (default 0):
%                               - 0 for the a zero matrix
%                               - 1 for A^t*y
%                               - give your own matrix initialization
%           -> itmaxIRL    : max iterations outer loop (default 10)
%           -> itmaxFista  : max iterations inner FISTA loop (default 10)
%           -> xtol        : tolerance on the relative difference btw
%                            to successive iterates (default 1e-5)
%           -> costFuncTest: boolean to activate the computation of the
%                            cost (default 0)
%           -> verbose     : boolean for verbose mode (default 0)
%
%
%  Ouputs - SRX : cell array of the Super-Resolved reconstructions
%                 obtained for each lambda
%         - infos : structure with the following fields:
%                    - lambda : vector of lambdas used
%                    - PSF    : PSF used
%                    - x0     : choosen initialisation
%                    - time   : vector of processing time for each recons
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

% === +++++++++++++++++++++++++++++++++++++++ === %
% === Set default values, warnings and errors === %
% === +++++++++++++++++++++++++++++++++++++++ === %
if size(Iacq,1) ~= size(Iacq,2), error('This algorithm was designed for square acquisitions only'); end
if isfield(options,'coefEch') && round(options.coefEch) ~= options.coefEch, error('The upsizing coefficient has to be an integer'); end
if ~isfield(options,'coefEch'), warning('The upsizing coefficient has not been indicated. Default value is 2'); options.coefEch = 2; end
if ~isfield(options, 'x0'),warning('The choice of X_0 has not been indicated. Default value is the zero matrix'); options.x0=0;end
if ~isfield(options,'itmaxIRL'),warning('The maximum number of iterations for the IRL1 loop has not been indicated. Default value is 10'); options.itmaxIRL=10; end
if ~isfield(options,'itmaxFista'), warning('The maximum number of iterations for the Fista loop has not been indicated. Default value is 100'); options.itmaxFista=100; end
if ~isfield(options,'xtol'), options.xtol=1e-5; end;
if ~isfield(options, 'costFuncTest'), options.costFuncTest = 0;end;
if ~isfield(options,'lambda'), error('Please indicate at least one lambda to be tested'); end
if ~isfield(options,'verbose'), options.verbose=0; end

SRSize = size(Iacq,1)*options.coefEch;
SRX = cell(1,size(options.lambda,2));
infos.lambda = options.lambda;

% === ++++++++++++++++++++++ === %
% === Convolution Parameters === %
% === ++++++++++++++++++++++ === %
PSF = options.PSF;
PSFfft = fft2(fftshift(PSF));
PSFconj = conj(PSFfft);
infos.PSF=PSF;

% === +++++++++++++++++ === %
% === Downsizing Matrix === %
% === +++++++++++++++++ === %
Mech = zeros(size(Iacq,1), SRSize);
for i = 0 : size(Iacq,1)-1
    Mech(i+1, 1+options.coefEch*i : options.coefEch+options.coefEch*i) = 1;
end
MechT = Mech';

% === +++++++++++++++ === %
% === Column A_i norm === %
% === +++++++++++++++ === %
if ~isfield(options,'normai')
    normai= ComputeNorm_ai(options.PSF,options.coefEch);
else
    normai=options.normai;
end

% === ++++++++++++++++++++ === %
% ===  Lipchitz Constant   === %
% === ++++++++++++++++++++ === %
lipCst = max(norm(Mech*MechT)*abs(PSFfft(:)))^2;
gamma = 1/lipCst;

% === ++++++++++++++++++ === %
% === X_0 Initialisation === %
% === ++++++++++++++++++ === %

if isscalar(options.x0)
    if options.x0 == 0
        xInit = zeros(SRSize);
        initText = 'a matrix of zeros';
    elseif  options.x0 == 1
        xInit = real(ifft2(conj(PSFfft).*fft2(MechT*Iacq*Mech)));
        initText = 'the retroprojected image';
    else
        error('Wrong choice for X_0 initialisation')
    end
else
    if size(options.x0,1) == SRSize && size(options.x0,2) == SRSize
        xInit = options.x0;
        initText = 'defined by user';
    else
        error('Wrong size for X_0 initialisation')
    end
end
infos.x0 = xInit;

% ===================== +++++++++++++++++++++++ =================== %
% ===================== Loop Algorithm IRL1-FBS =================== %
% ===================== +++++++++++++++++++++++ =================== %
if options.verbose == 1
    fprintf(' ________________________________________________\n');
    fprintf('|                                                |\n');
    fprintf('|                   SMLM-CEL0                    |\n');
    fprintf('|________________________________________________|\n\n');
    fprintf('Main parameters:                                 \n');
    fprintf('                 - %i iterations for outer loop \n',options.itmaxIRL);
    fprintf('                 - %i iterations for inner loop \n',options.itmaxFista);
    fprintf('                 - %i lambdas to be tested \n',size(options.lambda,2));
    fprintf('                 - Size of the acquisition %i x %i pixels \n',size(Iacq,1),size(Iacq,1));
    fprintf('                 - Size of the reconstruction %i x %i pixels \n', SRSize, SRSize);
    fprintf('                 - Initialisation of X is %s \n\n\n\n',initText);
end

% === +++++++++++++++++++++ === %
% === Loop over the lambdas === %
% === +++++++++++++++++++++ === %
for nbLam = 1:length(options.lambda)
    
    lambda =  options.lambda(nbLam);  % get the current tested lambda
    xIt = xInit;                      % initial point
    
    % === +++++++++ === %
    % === Loop IRL1 === %
    % === +++++++++ === %
    tic
    for it2 = 1:options.itmaxIRL
        xOld2=xIt;
        y = xIt;
        tk = 1;
        absx = abs(xIt);
        weights = (absx>=eps).*(absx<(sqrt(2*lambda)./normai)).*(-(normai.^2).*absx+(sqrt(2*lambda).*normai)) + ...
                  (absx<eps).*(sqrt(2*lambda).*normai);
        
        % === ++++++++++ === %
        % === Loop Fista === %
        % === ++++++++++ === %
        for it = 1: options.itmaxFista
            xOld = xIt;           
            told = tk;
            
            yfft = fft2(y);
            gradientF = real(ifft2(PSFconj.*fft2(MechT*(Mech * real(ifft2(PSFfft .* yfft)) * MechT - Iacq)*Mech)));
            xIt = y - gamma*gradientF;
            
            % === Prox L1 + positivity === %
            xIt = max((xIt-gamma*weights),0);
            tk = 0.5*(1+sqrt(1+4*tk^2));
            mu = (told-1)/tk;
            y = xIt + mu * (xIt - xOld);
            
            if options.costFuncTest == 1
                costOld =  0.5*norm(Mech*real(ifft2(PSFfft.*fft2(xOld)))*MechT-Iacq,2)^2 + norm(weights.*xIt, 1);
                costIt = 0.5*norm(Mech*real(ifft2(PSFfft.*fft2(xIt)))*MechT-Iacq,2)^2 + norm(weights.*xIt, 1);
                diffCout = abs(costIt - costOld)/abs(costIt+eps);
                if diffCout < options.xtol && norm(xIt(:)-xOld(:))/(norm(xOld(:))+eps) < options.xtol
                    break;
                end           
            elseif norm(xIt(:)-xOld(:))/(norm(xOld(:))+eps) < options.xtol,
                break;
            end
        end
        
        if options.costFuncTest == 1
            costOld = 0.5*norm(Mech*real(ifft2(PSFfft.*fft2(xOld2)))*MechT-Iacq,2)^2 +  sum(sum(lambda - 0.5 * normai.^2 .* (absx - sqrt(2*lambda)./normai).^2 .* (absx <= sqrt(2*lambda)./normai) ));
            costIt = 0.5*norm(Mech*real(ifft2(PSFfft.*fft2(xIt)))*MechT-Iacq,2)^2 +  sum(sum(lambda - 0.5 * normai.^2 .* (absx - sqrt(2*lambda)./normai).^2 .* (absx <= sqrt(2*lambda)./normai) ));
            diffCout = abs(costIt - costOld)/abs(costIt+eps);
            if diffCout < options.xtol && norm(xIt(:)-xOld2(:))/(norm(xOld2(:))+eps) < options.xtol
                break;
            end
        elseif norm(xIt(:)-xOld2(:))/(norm(xOld(:))+eps) < options.xtol,
            break;
        end
    end
    
    % === ++++++++++++++++++++ === %
    % === End of the algorithm === %
    % === ++++++++++++++++++++ === %
    SRX{nbLam} = xIt;
    infos.time(nbLam) =  toc;
    
    if options.verbose == 1
        fprintf(' _________________________________________\n');
        fprintf('                                         \n');
        fprintf('  Reconstruction ended for lambda = %i  \n', options.lambda(nbLam));
        fprintf('              Time = %.1f s              \n', infos.time(nbLam));
        fprintf('          %i non-zero pixels found      \n', sum(xIt(:)>eps));
        fprintf(' _________________________________________\n\n')
    end
end

if options.verbose == 1
    fprintf('+++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf('+                                         +\n');
    fprintf('+          End of the function            +\n');
    fprintf('+           Total time = %.1f s           +\n', sum(infos.time));
    fprintf('+                                         +\n');
    fprintf('+++++++++++++++++++++++++++++++++++++++++++\n')
end
end

