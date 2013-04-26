function r = synthesisNorm(mixture, mask, fRange, winLength, fs) 
% synthesize a waveform signal from a mixture signal and a binary mask.
% It amounts to an inverse filter of the gammatone filterbank. 
% The first two variables are required. 
% fRange: frequency range.
% winLength: window (frame) length in samples
% fs: sampling frequency.
% Written by ZZ Jin, adapted by DLW in Jan'07 and JF Woodruff in Nov'08

if nargin < 3
    fRange = [80, 5000]; % default frequency range in Hz
end
if nargin < 4
    winLength = 320;    % default window length in sample points which is 20 ms for 16 KHz sampling frequency
end
if nargin < 5
    fs = 16000;     % default sampling frequency
end

filterOrder = 4;     % filter order
winShift = winLength/2;     % frame shift rate (default is 1/2 frame length)
increment = winLength/winShift;     % special treatment for beginning frames
sigLength = length(mixture);
r = zeros(1,sigLength);
[numChan,numFrame] = size(mask);     % number of channels and time frames

for i = 1:winLength     % calculate a raised cosine window
    coswin(i) = (1 + cos(2*pi*(i-1)/winLength - pi))/2;
end

temp1 = gammatoneNorm(mixture, numChan, fRange, fs);          % first pass through gammatone filterbank

phase(1:numChan) = zeros(numChan,1);        % initial phases
erb_b = hz2erb(fRange);       % upper and lower bound of ERB
erb = [erb_b(1):diff(erb_b)/(numChan-1):erb_b(2)];     % ERB segment
cf = erb2hz(erb);       % center frequency array indexed by channel
b = 1.019*24.7*(4.37*cf/1000+1);       % rate of decay or bandwidth

% Generating gammatone impulse responses with middle-ear gain normalization
gL = 2048;      % gammatone filter length or 128 ms for 16 kHz sampling rate
gt = zeros(numChan,gL);
tmp_t = [1:gL]/fs;
for c = 1:numChan
    gain = 1/3*(2*pi*b(c)/fs).^4;    % loudness-based gain adjustments
    gt(c,:) = gain*fs^3*tmp_t.^(filterOrder-1).*exp(-2*pi*b(c)*tmp_t).*cos(2*pi*cf(c)*tmp_t+phase(c));    
end

% Gammatone response normalization
f=cf(3):cf(end-3); % the approximate plateu region
gtf = zeros(numChan,length(f));
for c=1:numChan
    gtf(c,:) = (1/3*(2*pi)^4)*factorial(filterOrder-1)*1/2*(b(c)./sqrt((2*pi*(f-cf(c))).^2+(2*pi*b(c))^2)).^filterOrder;
    gtf(c,:) = gtf(c,:).^2;
end
factor = sqrt(mean(sum(gtf)));
gt = gt/factor;

d = 2^ceil(log2(sigLength+gL));         % used in second path gammatone filtering
for c = 1:numChan    
    temp1(c,:) = fliplr(temp1(c,:));    % time reverse filter output
    temp2 = fftfilt(gt(c,:),temp1(c,:));    % second pass filtering via FFTFILT
    temp1(c,:) = fliplr(temp2(1:sigLength));    % time reverse again
    
    weight = zeros(1, sigLength);       % calculate weighting
    for m = 1:numFrame-increment/2+1      % mask value can be binary or rational           
        startpoint = (m-1)*winShift;
        if m <= increment/2                % shorter frame lengths for beginning frames or zero padding
            weight(1:startpoint+winLength/2) = weight(1:startpoint+winLength/2) + mask(c,m)*coswin(winLength/2-startpoint+1:end);
        else 
            weight(startpoint-winLength/2+1:startpoint+winLength/2) = weight(startpoint-winLength/2+1:startpoint+winLength/2) + mask(c,m)*coswin;
        end
    end
    
    r = r + temp1(c,:).*weight;
end