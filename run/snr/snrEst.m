function [masks, snr] = snrEst(sig,params)
if nargin < 2
    params = struct;
end
if isfield(params,'ebm')
    ebm = params.ebm;
end

if isfield(params,'tmpDir')
    tmpDir = params.tmpDir;
else
    tmpDir = pwd;
end

snr.readme = 'Estimated SNR values';
masks.readme = 'Estiamted IBMs';

wavnorm = 3.276797651191444e+004;
fs = 16000;
numChan = 128;
fRange = [50 8000];
winLength = 20*fs/1000;

  
if ~exist('ebm','var')
    fprintf(1,'\nPass a binary mask\n');
    return;
end

if ischar(sig)
    sig = round(wavread(sig)*wavnorm);
end

if ischar(ebm)
    ebm = load(ebm);
end

if isfield(params,'broadband')
    engy = cochleagram(gammatoneSNR(sig,numChan,fRange,fs),winLength);
    
    target = sum(sum(engy .* ebm));
    noise = sum(sum(engy)) - target;
    snr.filtered = 10*log10(target) - 10*log10(noise);
    
    ccgmEn = sum(sum(engy));
    diffEn = 2*sum(sig.^2) - ccgmEn;
    k=10^(snr.filtered/10);
    cleanEn = (k/(1+k))*ccgmEn;
    snr.broadband = 10*log10(cleanEn/(ccgmEn-cleanEn+max(0,diffEn)));
end
end
