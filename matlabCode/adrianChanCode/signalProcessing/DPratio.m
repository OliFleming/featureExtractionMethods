% DPratio This computes the maximum-to-minimum drop in power density ratio.
%
% DP = DPratio(f,p)
%
% Author Adrian Chan
%
% This computes the maximum-to-minimum drop in poewr density ratio (DPR).
% The DPR is the ratio between the highest mean power density and lowest
% mean power density between 35 and 600 Hz. The mean power density is
% computed by averaged 13 consecutive points in the EMG power spectrum. The
% DPR is used to indicate whether the power spectrum was adequately peaked
% in the expected frequency range and can be used to detect the absence of
% EMG activity. It also ensures that the EMG spectrum drops off for higher
% frequencies, which would enable detection of high frequency noise and
% aliasing due to undersampling.
%
% Reference: Sinderby C, Lindstrom L, Grassino AE, "Automatic assessment of
% electromyogram quality", Journal of Applied Physiology, vol. 79, no. 5,
% pp. 1803-1815, 1995.
%
% Inputs
%    f: frequencies (Hz)
%    p: power spectral density values
%
% Outputs
%    DP: maximum-to-minimum drop in power density (can be converted to dB 
%        by using 10*log10(SMratio(f,P))
%
% Modifications
% 09/09/21 AC First created.
function DP = DPratio(f,p)

debugmode = false;

% remove frequencies above 600 Hz
index_below_600 = find(f <= 600);
f = f(index_below_600);
p = p(index_below_600);

% average PSD over N points using N/2 points before and after
N = 13;
b = ones(N,1)/N;
a = 1;
mean_psd = filter(b,a,[p;zeros(floor(N/2),1)]);
mean_psd = mean_psd(floor(N/2) + (1:length(f)));

if debugmode == true
    figure
    plot(f,p,f,mean_psd);
    xlabel('f');
    ylabel('PSD');
    legend('PSD','Averaged PSD');
	title('DP ratio')
end

index_f_above_35 = find(f > 35);

highest_mean_psd = max(mean_psd(index_f_above_35));
lowest_mean_psd = min(mean_psd(index_f_above_35));

if debugmode == true
    f_highest_mean_psd = f(find(mean_psd == highest_mean_psd));
    f_lowest_mean_psd = f(find(mean_psd == lowest_mean_psd));
    hold on, plot(f_highest_mean_psd,highest_mean_psd,'ro'), hold off
    hold on, plot(f_lowest_mean_psd,lowest_mean_psd,'rs'), hold off
end

DP = highest_mean_psd/lowest_mean_psd;