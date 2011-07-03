function [y, errF] = forceFourierMagnitudeAndPhase(x, xSize, FMod, supportF, ...
        complexVarsFlag, phaseLo, phaseHi)
    % FORCEFOURIERMAGNITUDEANDPHASE - applies known Fourier magnitude and phase
    % interval (if provided) to x.



    % Copyright 2008-2011 Eli Osherovich.
    

    
    if nargin < 6
        phaseLo = [];
        phaseHi = [];
    end
    
    % To speed-up some computations we keep certain variables persisten in
    % the memory.
    persistent prevPhaseLo prevPhaseHi dirLo dirHi dirMid
    
    phaseLoCHANGE = ndims(prevPhaseLo) ~= ndims(phaseLo) || ...
        any(size(prevPhaseLo) ~= size(phaseLo)) || any(prevPhaseLo(:)~=phaseLo(:));
    phaseHiCHANGE = ndims(prevPhaseHi) ~= ndims(phaseHi) || ...
        any(size(prevPhaseHi) ~= size(phaseHi)) || any(prevPhaseHi(:)~=phaseHi(:));
    
    % Direction of the half-plane "above" the lower phase bound.
    if phaseLoCHANGE
        dirLo = exp(1i*(phaseLo+pi/2));
        prevPhaseLo = phaseLo;
    end
    
    % Direction of the half-plane "below" the upper phase bound.
    if phaseHiCHANGE
        dirHi = exp(1i*(phaseHi-pi/2));
        prevPhaseHi = phaseHi;
    end
    
    % Direction of the middle half-plane (towards phaseHi).
    % This half-plane bisects the angle between phaseLo and phaseHi hence,
    % any copmlex number lying in the "upper" (positive) half is close to
    % phaseHI. Consequently, any complex number lying the "lower"
    % (negative) half is closer to phaseLo.
    if phaseLoCHANGE || phaseHiCHANGE
        dirMid = exp(1i*((phaseLo+phaseHi)/2 + pi/2));
    end
    
    
    % Fourier transform operator (with padding).
    F = UnitaryDFT_eo(xSize);
    X = F*x;
    
    % In general, Fouier domain data is available only for some
    % frequencies (defined by supportF).
    XS = X(supportF);
    XSPhase = angle(XS);
    
    
    % Enforce phase constraints.
    if ~(isempty(dirLo) || isempty(dirHi) || isempty(dirMid))
        % Projections on the half-planes defined by dirLo, dirHi, and dirMid.
        dirLoProj = real(conj(dirLo).*XS);
        dirHiProj = real(conj(dirHi).*XS);
        dirMidProj = real(conj(dirMid).*XS);
        
        % Violations of the upper bound (phaseHi).
        violHiIdx = (dirMidProj > 0) & (dirHiProj < 0);
        
        % Violations of the lower bound (phaseLo).
        violLoIdx = (dirMidProj <=0) & (dirLoProj < 0);
        
        % Force phase bounds.
        XSPhase(violHiIdx) = phaseHi(violHiIdx);
        XSPhase(violLoIdx) = phaseLo(violLoIdx);
    end
    
    % New XFS with correct Fourier magnitued and phase withing the bounds.
    XSF = FMod.*exp(1i*XSPhase);
    X(supportF) = XSF;
    
    % Error in the Fourier domain.
    errF = norm(XSF - XS)^2;
    
    % Convert back to the object domain. We use the fact that F is
    % orthogonal, i.e., inv(F) = F'.
    if complexVarsFlag
        y = F'*X;
    else
        y = real(F'*X);
    end
