function [sig_noisy, m_photons] = generateNoisySig_Poisson_SNR(squared_sig, target_snr)
% generate a noisy signal with Poisson noise of given SNR



% Copyright 2008-2010 Eli Osherovich.


% fast exit if INF snr requested (it corresponds to uncontaminated signal)
if Inf == target_snr
    sig_noisy = sqrt(squared_sig);
    m_photons = Inf;
    return;
end


% snr tolerance
snr_tol = 0.1;

% accepted snr range
min_snr = target_snr - snr_tol;
max_snr = target_snr + snr_tol;


%% run bisection method to obtain the noise of requested SNR

% find an appropriate interval
[a_photons, b_photons] = find_interval(squared_sig, target_snr);

% run bisection method
while b_photons > a_photons
    m_photons = round((a_photons + b_photons)/2);
    [snr_curr, sig_noisy] = calc_noisy_sig(squared_sig, m_photons);
    if snr_curr > max_snr
        b_photons = m_photons;
    elseif snr_curr < min_snr
        a_photons = m_photons;
    else
        return;
    end
end




function [a_photons, b_photons] = find_interval(squared_sig, target_snr)

maxIter = 100;


% start with guessing a good photon number
a_photons = 1e9;
b_photons = 1e12;

snr_a = calc_noisy_sig(squared_sig, a_photons);
snr_b = calc_noisy_sig(squared_sig, b_photons);

if target_snr < snr_a
    for i = 1:maxIter
        snr_new = calc_noisy_sig(squared_sig, a_photons/2);
        if target_snr > snr_new
            b_photons = a_photons;
            a_photons = a_photons/2;
            return
        else
            a_photons = a_photons/2;
        end
    end
    if i == maxIter
        warning('generateNoisySig_Poisson_SNR:maxIter', 'Could not find valid interval'); 
    end
elseif target_snr > snr_b
    for i = 1:maxIter
        snr_new = calc_noisy_sig(squared_sig, b_photons*2);
        if target_snr < snr_new
            a_photons = b_photons;
            b_photons = 2*b_photons;
            return
        else
            b_photons = 2*b_photons;
        end
    end
else
    return
end

    

function [snr, sig_Noisy] = calc_noisy_sig(squared_sig, n_photons)

% keep current rand stream
old_rstream = RandStream.getDefaultStream;

% use a separate RandStream
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', 0));

% convert grayscale into photons
a = n_photons/sum(squared_sig(:));
squared_sig_photons = squared_sig*a;

squared_sig_photons_Noisy = poissrnd(squared_sig_photons);
% violations
violations = squared_sig_photons_Noisy < 0;
squared_sig_photons_Noisy(violations) = 0;

squared_sig_Noisy = squared_sig_photons_Noisy/a;

sig = sqrt(squared_sig);
sig_Noisy = sqrt(squared_sig_Noisy);
noise = sig_Noisy - sig;

snr = 20*log10(norm(sig(:))/norm(noise(:)));

% restore original default RandStream
RandStream.setDefaultStream(old_rstream);
