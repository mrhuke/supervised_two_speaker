function [snr gain]= synthesis_SNR(msk,ibm,sig,fRange,winLen,fs)

addpath(genpath('../'));

if nargin<=3
    fRange = [50 8000];
    winLen = 320;
    fs = 16000;
end

estSig = synthesis(sig, msk, fRange, winLen, fs);
iSig = synthesis(sig, ibm, fRange, winLen, fs);
allOneMask = ibm;
allOneMask(:,sum(ibm)>0) = 1;
initSig = synthesis(sig, allOneMask, fRange, winLen, fs);

snr = 10*log10(iSig.^2/(estSig-iSig).^2);
initSNR = 10*log10(iSig.^2/(initSig-iSig).^2);
gain = snr-initSNR;

