%Q1_Apart
length = 1;
frequency = 1000;
sigma = 1;
T = 0.1;
plot(digitalWave(length,frequency,sigma,T))
xlabel('time(s)')
ylabel('x(t)')
title('Random Digital Wave')

%Q1_B
length = 1;
frequency = 1000;
sigma = 1;
T = 0.1;
numberOfSignals = 1000;
[meanCorr,] = autocorrelation(length,frequency,sigma,T,numberOfSignals);

plot(meanCorr)
xlabel('signals!!!(There are 1000 signals)')
ylabel('meanCorrelation')
title('part B')

%Q1_C
% sweep from sigma = 0.1 to 1 Dont run this at all :)))) just read the
% results
length = 1;
frequency = 1000;
numberOfSignals = 1000;
T = 0.1;
for i = 1:5
    sigma = i/5;
    [meanCorr,] = autocorrelation(length,frequency,sigma,T,numberOfSignals);
    figure
    plot(meanCorr)
    xlabel('signals!!!(There are 1000 signals)')
    ylabel('R_xx(t)');
    title(['sigma = ',num2str(sigma)]);
end

for i = 1:5
    T = i/10;
    figure
    [meanCorr,] = autocorrelation(length,frequency,sigma,T,numberOfSignals)
    figure
    plot(meanCorr)
    xlabel('signals!!!(There are 1000 signals)')
    ylabel('R_xx(t)');
    title(['T = ',num2str(i/10)]);
end

%Q1_D
length = 1;
frequency = 1000;
numberOfSignals = 1000;
T = 0.1;
sigma = 1;
[meanCorr,allsignals,correlations] = autocorrelation(length,frequency,sigma,T,numberOfSignals);
energy_viaRxx = abs(sdfViaCorrSignal(meanCorr));
plot(energy_viaRxx)
xlabel('frequency')
ylabel('Standard energy spectrum calculted via R_x_x')
title('part D')

%Q1_E
% [meanCorr,allsignals,correlations] = autocorrelation(length,frequency,sigma,T,numberOfSignals);
energy = zeros(size(allsignals,1),size(allsignals,2));
for i=1:size(allsignals,1)
    extractedSignal = allsignals(i,:);
    fourierTrans = fftshift(fft(extractedSignal,size(extractedSignal,2)));
    energy(i,:) = fourierTrans.*conj(fourierTrans);
end
Energy_fourierTrans = mean(energy);
figure
plot(abs(Energy_fourierTrans));
xlabel('frequency')
ylabel('Standard energy spectrum calculted via X(f)')
title('part E')
MSE = abs(energy_viaRxx - Energy_fourierTrans).^2;
MSE_mean = mean(MSE);
plot(MSE)
xlabel('frequency')
ylabel('MSE of (Sdf calculated from R_x_x & Sdf calculted from X(f)')
title('part E')


%Q1_F
length = 1;
frequency = 1000;
T = 0.1;
sigma = 1;
numberOfSignals = 1000;

eachSignalSize = length*frequency;  
allSignals = zeros(numberOfSignals,eachSignalSize); % we can change it by fixing of the outputSize!
for i = 1:numberOfSignals
    allSignals(i,:) = digitalWave(length,frequency,sigma,T);
end

sinDomain =  0 : 1/frequency : (size(allSignals,2)-1)/frequency; % syncing both sampling frequency & signalDomainNumber
sinus = sin(2*pi*200*sinDomain);
out = dotDigitalWave(sinus,allSignals);
figure
plot(abs(out))
xlabel('frequency')
ylabel('Amplitude(|H(f)|)')
title('part F')
figure
plot(abs(sdfViaCorrSignal(out)))
xlabel('frequency')
ylabel('Standard energy spectrum')
title('part F')


%Q2_a
[audioData,fs] =  audioread('clip.wav');
plot(audioData);
xlabel('n');
ylabel('AudioFile[n]');
title('AudioFileProcessing');
[audioDataSdf,audioDataFourierTransform] = sdf(transpose(audioData));
figure
plot(abs(audioDataFourierTransform))
xlabel('frequency');
ylabel('AudioDataAfterFilterFourierTransform[n]');
title('AudioFileProcessing');
figure
plot(abs(audioDataSdf))
xlabel('frequency');
ylabel('AudioDataAfterFilterEnergy[n]');
title('AudioFileProcessing');

%Q2_b

% first verifying linearinty of the channel
domain = 0:0.01:2*pi;
inp1 = sin(domain);
inp2 = cos(domain);
totSignal = inp1+inp2;

out1 = Channel(inp1);
out2 = Channel(inp2);
totOut = Channel(totSignal);

figure
subplot(2,1,1)
plot(domain,inp1)
xlabel('n')
ylabel('inp1[n]')
title('firstInput')
xlim([0 2*pi])
grid on
subplot(2,1,2)
plot(domain,inp2)
xlabel('n')
title('secondInput')
ylabel('inp2[n]')
xlim([0 2*pi])
grid on
figure
subplot(3,1,1)
plot(domain,out1+out2,'b*')
xlabel('n')
ylabel('Channel(inp1) + Channel(inp2)')
title('LinearityCheck')
grid on
xlim([0 2*pi])
subplot(3,1,2)
plot(domain,totOut,'r*')
title('LinearityCheck')
xlabel('n')
ylabel('Channel(inp1+inp2)')
xlim([0 2*pi])
grid on
subplot(3,1,3)
MSE_linearity = abs(totOut-(out1+out2)).^2;
plot(domain,MSE_linearity,'g*')
xlabel('n')
ylabel('MSE error[n]')
title('MSE error')
xlim([0 2*pi])
grid on

% second verifying timeVariant Check
shiftValue = 50;
x = 0:0.01:2*pi;
xShifted = [zeros(1,shiftValue) x];
y = Channel(x);
yShifted = Channel(xShifted);
figure
subplot(3,1,1)
plot(y,'r','linewidth',2)
title('The output of x')
ylabel('Channel(x)[n]')
xlabel('[n]')
subplot(3,1,2)
plot(yShifted,'r','linewidth',2)
title('The output of the shifted x')
xlabel('[n]')
ylabel('Channel(shiftedX)[n]')
MSE_timeInvariant = abs(yShifted(shiftValue + 1:end)-y).^2;
subplot(3,1,3)
plot(MSE_timeInvariant)
title('MSE of Channel(x) & Channel(shiftedX)')
xlabel('n')
ylabel('MSE of shifted & MainData')
%Q2_c groupDelay & phaseDelay! 
% giving System delta input with the size of given clip.wmv
delta = [1 zeros(1,size(transpose(audioData),2)-1)];
h = Channel(delta);
n = 0:size(h,2)-1;
% using duality in the time & frequency Domain.
% derrivation in frequency == multipication in time domain
groupDelay = real(fourierDiscreteStandardTransform(n.*h)./fourierDiscreteStandardTransform(h))/2*pi;

%syncing input among -pi til pi.
frequency = 22050; % Hz
f = (0:size(groupDelay,2)-1)*(frequency)/(size(groupDelay,2)-1) - frequency/2;
figure
plot(f,groupDelay)
xlabel('Frequency(Hz)')
ylabel('gd(f)')
title('groupDelay')
H = fourierDiscreteStandardTransform(h);
phase_delay = -angle(H)./(2*pi*f);
figure
plot(f,phase_delay)
xlabel('Frequency(Hz)')
ylabel('pd(f)')
title('Phase delay of channel')

%Q2_d
audioChannelOutput = Channel(audioData);
figure
plot(audioChannelOutput);
xlabel('n');
ylabel('Channel(AudioData)[n]')
title('AudioFileProcessing')
[audioChannelOutputSdf,audioChannelOutputFourierTransform] = sdf(transpose(audioChannelOutput));
figure
plot(abs(audioChannelOutputFourierTransform))
xlabel('frequency')
ylabel('audioChannelOutputFourierTransform')
title('FourierTransformAfterChannel')
figure
plot(abs(audioChannelOutputSdf))
xlabel('frequency')
ylabel('audioChannelAfterFilterEnergy')
title('Sdf after Channel')

%Q2_e
% use band pass filter for those who are not too little.
threshHoldValue = 0.0000001;
gonnaRemoveIndexes = find(abs(H(floor(size(H,2)/2):end)) < threshHoldValue*max(abs(H)));
% remove from the first one.
BoundeStop = gonnaRemoveIndexes(1);
% Equalizer Setting.
Equalizer = zeros(1,size(H,2));
Pass_band = floor(size(H,2)/2)- BoundeStop+1 : floor(size(H,2)/2)+ BoundeStop+1;
Equalizer(Pass_band) = H(Pass_band).^-1;
plot(f,abs(Equalizer))
xlabel('f')
ylabel('|H(f)|')
title('|H_e_q_u_a_l_i_z_e_r(f)|')
figure
plot(f,angle(Equalizer))
xlabel('f')
ylabel('phase(rad)')
title('angle(H_e_q_u_a_l_i_z_e_r(f))')

% checking the accurately of Equalizer
[audioData,fs] =  audioread('clip.wav');
audioDataTransform = fourierDiscreteStandardTransform(transpose(audioData));
audioChannelOutput = Channel(audioData);
afterEqualizer = fourierDiscreteStandardTransform(transpose(audioChannelOutput)).*Equalizer;
MSE = abs((afterEqualizer - audioDataTransform).^2);
plot(MSE)
xlabel('frequency')
ylabel('MSE')
title('Checking Accuracy of Equalizer')

%Q2_f
% input output for Equalizer
% I told all the steps done this section in the report file.
systemIdentification

%Q3_a,Q3_b
frequency = 22050;
length = 1000; % youConsider.
[hilbertSignal,zeroIndex,domain] = hilbert(length,frequency);
% removing Zero values in order to enhance the curve.
% contact Me!
stem(domain,hilbertSignal);
xlim([-0.001 0.001])
xlabel('n')
title('Hilbert Implulse Response')
ylabel('ImpluseResponse[n]')

%Q3_c
M = 500;
numberOfDatas = M;
cutDot = floor(numberOfDatas/2);
afterFilter = zeros(1,2*cutDot+1);
afterFilter(1:end) = hilbertSignal(1,zeroIndex-cutDot:zeroIndex+cutDot);
myDomain = domain(zeroIndex-cutDot:zeroIndex + cutDot);
shiftSize = floor(M/2) + 1;
shiftedHilbert = circshift(afterFilter,shiftSize);
shiftedDomain = myDomain + shiftSize/frequency;
% critical condition happened Data Around
hilbertTransform = fourierDiscreteStandardTransform(shiftedHilbert);
w = (0:size(shiftedHilbert,2)-1)*2*pi/(size(shiftedHilbert,2)-1) - pi;
hilbertTransform = hilbertTransform.*exp(1i*-0.5*w);
magnitude = abs(hilbertTransform);
angle1 = angle(hilbertTransform);
figure
subplot(2,1,1)
plot(magnitude)
xlabel('frequency');
ylabel('Amplitude');
xlim([-100 600])
title(['M = ',num2str(M)]);
subplot(2,1,2)
plot(angle1)
xlabel('frequency')
xlim([-100 600])
ylabel('Angle')

%Q3_d
M = 500;
numberOfDatas = M;
cutDot = floor(numberOfDatas/2);
afterFilter = zeros(1,2*cutDot+1);
afterFilter(1:end) = hilbertSignal(1,zeroIndex-cutDot:zeroIndex+cutDot);
myDomain = domain(zeroIndex-cutDot:zeroIndex + cutDot);
shiftSize = floor(M/2) + 1;
shiftedHilbert = circshift(afterFilter,shiftSize);
shiftedDomain = myDomain + shiftSize/frequency;
afterHammingFilter = shiftedHilbert.*hamming(domain(zeroIndex-cutDot:zeroIndex+cutDot),M,shiftSize);
myDomain = domain(zeroIndex-cutDot:zeroIndex + cutDot);
% shiftedHilbert = zeros(1,size(afterRect,2));
shiftedHilbert = myDomain + M/2;
figure

afterFilterTransform = fourierDiscreteStandardTransform(afterHammingFilter);
w = (0:size(afterFilterTransform,2)-1)*2*pi/(size(afterFilterTransform,2)-1) - pi;
afterFilterTransform = afterFilterTransform.*exp(1i*-0.5*w);

figure
subplot(2,1,1)
plot(abs(afterFilterTransform))
xlabel('frequency');
ylabel('Amplitude');
xlim([-100 600])
title(['Hamming M = ',num2str(M)]);
subplot(2,1,2)
plot(angle(afterFilterTransform))
xlabel('frequency');
xlim([-100 600])
ylabel('Angle');

% Last 
%% Bonus Part (Kaiser Window)
M = 300;
numberOfDatas = M;
cutDot = floor(numberOfDatas/2);
afterFilter = zeros(1,2*cutDot+1);
afterFilter(1:end) = hilbertSignal(1,zeroIndex-cutDot:zeroIndex+cutDot);
myDomain = domain(zeroIndex-cutDot:zeroIndex + cutDot);
shiftSize = floor(M/2) + 1;
shiftedHilbert = circshift(afterFilter,shiftSize);
shiftedDomain = myDomain + shiftSize/frequency;
kaiserVector = kaiser(M+1,7.5);
wvtool(kaiserVector)
kaiserShift = circshift(kaiserVector,shiftSize);
kaiserFilterTimeDomain = shiftedHilbert .* transpose(kaiserShift);
kaiserFilterFrequencyDomain = fourierDiscreteStandardTransform(kaiserFilterTimeDomain);
w = (0:size(kaiserFilterFrequencyDomain,2)-1)*2*pi/(size(kaiserFilterFrequencyDomain,2)-1) - pi;
kaiserFilterFrequencyDomain = kaiserFilterFrequencyDomain.*exp(1i*-0.5*w);

figure
plot(abs(kaiserFilterFrequencyDomain));
xlabel('frequency');
ylabel('Amplitude');
title('Frequency Answer(Kaiser Filter)');
xlim([-100 400])
figure
plot(angle(kaiserFilterFrequencyDomain));
xlabel('frequency');
ylabel('Angle');
xlim([-100 400])
title('Frequency Answer(kaiser Filter)');

function out = dotDigitalWave(digitalWave,signals)
signalDotDigitalWave = digitalWave .* signals;
[numberOfSignals,eachSignalSize] = size(signals);
correlation = zeros(numberOfSignals,eachSignalSize);
for i = 1:numberOfSignals
    for j = 1:eachSignalSize
    correlation(i,j) = Correlation(signalDotDigitalWave(i,:),circshift(signalDotDigitalWave(i,:),j-1));
    end   
end
mean_correlation = mean(correlation,1);
out = mean_correlation;
end

function out = sdfViaCorrSignal(signal)
out = fftshift(fft(signal,size(signal,2)));
end

function [out,allSignals,correlation] = autocorrelation(length,frequency,sigma,T,numberOfSignals)
eachSignalSize = length*frequency;  % length & 
allSignals = zeros(numberOfSignals,eachSignalSize); % we can change it by fixing  the outputSize!
for i = 1:numberOfSignals
    allSignals(i,:) = digitalWave(length,frequency,sigma,T);
end
correlation = zeros(numberOfSignals,eachSignalSize);
for i = 1:numberOfSignals
    for j = 1:eachSignalSize
    correlation(i,j) = Correlation(allSignals(i,:),circshift(allSignals(i,:),j-1));
    end   
end
out = mean(correlation,1);
end

function out = digitalWave(length,frequencySampling,sigma,T)
Ddomain = rand * T;
zeroPartNumbers = floor(Ddomain*frequencySampling);
oneIntervalDomainSamples = floor(T*frequencySampling);
discreteDomain = floor(length*frequencySampling);
numberOfIntervals = floor((discreteDomain - zeroPartNumbers)/oneIntervalDomainSamples);
discreteOutputSignal = zeros(1,discreteDomain);
discreteOutputSignal(1,1:zeroPartNumbers) = 0; % during zero part it is zero.
startIndex = zeroPartNumbers + 1;
for i = 1:numberOfIntervals
    a_i = normrnd(0,sigma);
    if(i~=numberOfIntervals)
    discreteOutputSignal(1,startIndex:startIndex + oneIntervalDomainSamples) = a_i;
    startIndex = startIndex + oneIntervalDomainSamples + 1;
    else
    startIndex = startIndex + 1;
    discreteOutputSignal(1,startIndex:end) = a_i;
    end
end
out = discreteOutputSignal;
end

function output = Correlation(firstSignal,secondSignal)
output = sum(firstSignal.*conj(secondSignal));
end


function out = fourierDiscreteStandardTransform(signal)
out = fftshift(fft(signal,size(signal,2)));
end

function out = hamming(input,M,shift)
temp = 0.54-0.46*cos(2*pi*input/M);
out = circshift(temp,shift);
end

function [out,zeroIndex,domain] = hilbert(length,frequency)
domain = (-1)*length/2 + 1/(2*frequency) :1/frequency:length/2 + 1/(2*frequency);
signal = zeros(1,size(domain,2));
zeroInit = ceil(size(domain,2)/2);
signal(1:end) = 1./(pi*domain(1:end));
zeroIndex = zeroInit;
out = signal;
end

function [out,Tsignal] = sdf(signal)
Tsignal = fourierDiscreteStandardTransform(signal);
out = Tsignal.*conj(Tsignal);
end
