% phase one.
samplePicture = imread('1.gif');
imshow(samplePicture);
title('original:512 & 512');
shrinkRate = 0.125;
beforeCoding = imresize(samplePicture,shrinkRate)
[Iindex , Jindex] = size(beforeCoding);
imshow(imresize(beforeCoding,1/shrinkRate));
title('Image:beforeCoding(64&64 rescaled to 512&512)');
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code)
sendingBinaryCode = createStringBinary(codedImage,0)
% decoding the image
afterDecoding = extract(sendingBinaryCode,MCode2Label,Iindex,Jindex,0)
modifieSizeImage = imresize(afterDecoding,1/shrinkRate);
imshow(modifieSizeImage);
title('Image:AfterDecoding(64&64 rescaled to 512&512)');

% phase two -> A
samplePicture = imread('1.gif');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
imshow(imresize(samplePictureResized,1/shrinkRate));
title('Image(resized 64&64 to 512&512):before module/demodulation');
[Iindex , Jindex] = size(samplePictureResized);
[labels,ps] = compile(samplePictureResized); % ps is sorted from more to less
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,0));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#'); % for noisy sending
% modulating 
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
% demoduling 
binaryDemodulatedSignal = demodule(analModulatedSignal,0,fs,fc,Ts,flags,0,0);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,0);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
imshow(modifieSizeImage);
title('Image(resized 64&64 to 512&512):after module/demodulation')

% phase two -> B
samplePicture = imread('1.gif');
%imshow(samplePicture);
%title('512 & 512');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
imshow(samplePictureResized);
%title('64 & 64');
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
%entering channel
%demodInput = Channel(analModulatedSignal,fc,fs,W)
% demoduling 
binaryDemodulatedSignal = demodule(analModulatedSignal,1,fs,fc,Ts,flags,1);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
imshow(modifieSizeImage);

% phase two -> C
samplePicture = imread('1.gif');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
[signalFourierTransform,finalDomain] = standardFourierTransform(analModulatedSignal,fs);
figure()
plot(finalDomain,abs(signalFourierTransform));
title('fourier transform of modulatedSignal');
xlabel('frequency');
ylabel('magnitude');
[audioDataSdf,~] = sdf(analModulatedSignal,fs);
tempCenter = find(finalDomain >= 10000);
centerFreq = tempCenter(1);
zeroTemp = find(finalDomain >= 0);
zero = zeroTemp(1);
powerAbs = abs(audioDataSdf);
totalEnergy = sum(powerAbs(zero:end));
power99 = totalEnergy*0.99;
freqIterator = 100;% 1kHz
boundedEnergy = calEnergy(powerAbs,centerFreq,freqIterator);
while boundedEnergy < power99
    freqIterator = freqIterator + 1000;
    boundedEnergy = calEnergy(powerAbs,centerFreq,freqIterator);
end
startBound = centerFreq - freqIterator;
endBound = centerFreq + freqIterator;
normalizedStartBound = finalDomain(startBound:endBound);
BW99percent = (vpa(normalizedStartBound(end)) - vpa(normalizedStartBound(1))) % it is bandPass there is no need to /2;
cutSignal = zeros(1,size(normalizedStartBound,2));
cutSignal(1:end) = powerAbs(startBound:endBound);
figure
plot(finalDomain,powerAbs)
xlabel('frequency');
ylabel('DataEnergy(f)');
title('modulatedSignal');
figure
plot(normalizedStartBound,cutSignal)
xlabel('frequency');
ylabel('DataEnergy(f)');
title('AfterChannel');
hold on
plot(normalizedStartBound(end),0,'*r')
plot(normalizedStartBound(1),0,'*r')
figure
plot(finalDomain(zero:end),powerAbs(zero:end))
xlim([9000 11000])
xlabel('frequency');
ylabel('DataEnergy(f)');
title('EnergySpectrum');
hold on
plot(normalizedStartBound(end),0,'*r')
plot(normalizedStartBound(1),0,'*r')
finalDomain(zero);
finalDomain(end);
plot(finalDomain(end),0,'*g')
plot(finalDomain(zero),0,'*g')
xlim([0 51000])
%entering channel
demodInput = Channel(analModulatedSignal,fc,fs,W)
[signalFourierTransform,finalDomain] = standardFourierTransform(demodInput,fs);
figure()
plot(finalDomain,abs(signalFourierTransform));
title('fourier transform of demodulated signal')
xlabel('frequency');
ylabel('magnitude');
% demoduling 
binaryDemodulatedSignal = demodule(demodInput,1,fs,fc,Ts,flags,0,0);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
imshow(modifieSizeImage);
title('Image:(resized 64&64 to 512&512)-> demodulatedImage')

% phase two -> C -> calculating Average BW of all Signals
% for in pictures
numberofPictures = 40;
BWs = zeros(1,numberofPictures);
name = '1.gif';
for i = 1:numberofPictures
number = i;
dotIndex = find(name == '.');
name = char(string(num2str(number)) + "." + string(name(dotIndex+1:end)))
samplePicture = imread(name);
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,0));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
BWs(i) = BW(analModulatedSignal,fs);
end
vpa(mean(BWs))

% phase two -> D
samplePicture = imread('1.gif');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
imshow(samplePictureResized);
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
% entering channel
noiseSNR = [10^(-1) 10^(-50) 10^(-100) 10^(-150)];
errorP = zeros(1,4);
for i =1:4
noisyAnal = awgn(analModulatedSignal/150,noiseSNR(i))

% mysigPower = sum(analModulatedSignal.*conj(analModulatedSignal));
% noies

demodInput = Channel(noisyAnal,fc,fs,W)
% demoduling 
binaryDemodulatedSignal = demodule(demodInput*150,1,fs,fc,Ts,flags,0,0);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
figure
imshow(modifieSizeImage);
end

% phase two -> E
% for in pictures
numberofPictures = 40;
snrsbefore = zeros(1,numberofPictures);
snrafter = zeros(1,numberofPictures);
name = '1.gif';

for i = 1:numberofPictures
number = i;
dotIndex = find(name == '.');
name = char(string(num2str(number)) + "." + string(name(dotIndex+1:end)))
samplePicture = imread(name);
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
% entering channel
noisyAnal = awgn(analModulatedSignal/150,0.00001);
demodInput = Channel(150*noisyAnal,fc,fs,W);
% demoduling 
% snr before
snrsbefore(i) = snr(demodInput, demodInput - analModulatedSignal);
binaryDemodulatedSignal = demodule(demodInput,1,fs,fc,Ts,flags,0,0);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
% snr after
mysize = size(samplePicture);
snrafter(i) = snr(double(reshape(samplePicture,[1 mysize(1,1)*mysize(1,2)])),double(reshape(modifieSizeImage - samplePicture,[1 mysize(1,1)*mysize(1,2)])));

if( 10 >= i && i >= 6 )
    figure
    imshow(imresize(samplePictureResized,1/shrinkRate));
    title('mainImage');
    figure
    imshow(modifieSizeImage)
    title('noisyImage');
end

end
meanSnrBefore = mean(snrsbefore)
meanSnrAfter = mean(snrafter)

% phase two -> F
samplePicture = imread('1.gif');
%imshow(samplePicture);
%title('512 & 512');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
%imshow(samplePictureResized);
%title('64 & 64');
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
% entering channel

% todo
% todo
noiseVariance = [0.01 0.0001 0.000001 0.00000001 0.0000000001 0.000000000001];

for i = 1:5
% increasing Variance ...
noisyAnal = awgn(analModulatedSignal/(i*50),noiseVariance(3));
demodInput = Channel(noisyAnal,fc,fs,W);
% demoduling 
binaryDemodulatedSignal = demodule(demodInput*(i*50),1,fs,fc,Ts,flags,0,1);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
figure
imshow(modifieSizeImage);
title("noiseVarianceEffect");
end
%

% phase two -> G
samplePicture = imread('1.gif');
%imshow(samplePicture);
%title('512 & 512');
shrinkRate = 0.125;
samplePictureResized = imresize(samplePicture,shrinkRate);
[Iindex , Jindex] = size(samplePictureResized);
%imshow(samplePictureResized);
%title('64 & 64');
[labels,ps] = compile(samplePictureResized); % ps is sorted from less to more
codes = encoding(ps);
keySet = labels;
valueSet = codes;
Mlabel2Code = containers.Map(keySet,valueSet);
MCode2Label = containers.Map(valueSet,keySet); % coding the image
codedImage = code(samplePictureResized,Mlabel2Code);
sendingBinaryCode = char(createStringBinary(codedImage,1));

% characteristics
fs = 100000; % 100khz
Ts = 0.01; % 10ms
fc = 10000; % 10khz
W = 5000; % 5khz
%

flags = find(char(sendingBinaryCode) == '#');
analModulatedSignal = modulate(sendingBinaryCode,fs,fc,Ts);
% entering channel
noisyAnal = awgn(analModulatedSignal/70,0.001);
demodInput = Channel(noisyAnal,fc,fs,W);
% demoduling 
angle = [0 pi/12 pi/6 pi/3];
for i = 1:4
binaryDemodulatedSignal = demodule(demodInput*70,1,fs,fc,Ts,flags,angle(i),1);
extractedImage = extract(binaryDemodulatedSignal,MCode2Label,Iindex,Jindex,1);
modifieSizeImage = imresize(extractedImage,1/shrinkRate);
figure
imshow(modifieSizeImage);
title("phase = " + string((angle(i)/pi)*180));
end
% 

function [labels,ps] = compile(image) % get 64 & 64 image
[horiSize , veriSize] = size(image);
imageVector = reshape(image,[1 horiSize*veriSize]);
imageVector = sort(imageVector);

pixels = removeRepeat(imageVector);
Ps = zeros(1,length(pixels));

mysize = length(pixels);
for i = 1:mysize
    Ps(i) = repeatition(pixels(i),image)/(veriSize*horiSize);
end

keySet = pixels;
valueSet = Ps;
M = containers.Map(keySet,valueSet);

sortedPs = varon(sort(Ps));
sortedPixels = zeros(1,length(sortedPs));

convertedType = zeros(1,length(sortedPs));

for i = 1:length(pixels)
        convertedType(i) = pixels(i);
end

for i = 1:length(sortedPs)
    for j = 1:length(pixels)
            if((convertedType(j) ~= -1) && (M(convertedType(j)) == sortedPs(i)))
                sortedPixels(i) = convertedType(j);
                convertedType(j) = -1;
                break;
            end
    end
end

labels = sortedPixels;
ps = sortedPs;

end


function removed = removeRepeat(vector)
mySize = length(vector);
temp = vector(1,1);
result = vector(1,1);
for i = 2:mySize
    if(vector(i) ~= temp)
        result = [result vector(i)];
        temp = vector(i);
    end
end
removed = result;
end


function repeat = repeatition(number,mat)
[horiSize , veriSize] = size(mat);
counter = 0;
for i = 1:horiSize
    for j = 1:veriSize
        if(mat(i,j) == number)
            counter = counter + 1;
        end
    end
end
repeat = counter;
end


function data = encoding(p)
words = string(zeros(1,length(p)));
for i = 1:length(words)
    words(i) = "";
end

% sorting p
% you should use it before the main call
%p = sort(p);
%sortedP = zeros(1,length(p));
%for i = 0:length(p)-1
%   sortedP(1,i+1) = p(1,length(p)-i);
%end
%p = sortedP;
%

data = shannon(1,length(p),p,words);
end

function data =  shannon(firstIndex,endIndex,p,words)

if(firstIndex ~= endIndex)
    
for i=firstIndex:endIndex
    
    sum1 = sum(p(firstIndex:i));
    sum2 = sum(p(i+1:endIndex));
    
    if( sum1>=sum2 )
        
        sum3 = sum(p(firstIndex:i-1));
        sum4 = sum(p(i:endIndex));
        
        diff1 = abs(sum1 - sum2);
        diff2 = abs(sum3 - sum4);
        
        if(diff1>diff2)
            index = i-1;
        else
            index = i;
        end
        
        break;
    end
end

words(firstIndex:index) = words(firstIndex:index) + "0"; % * is zero
words(index+1:endIndex) =  words(index+1:endIndex) + "1"; % # is one

if(index > firstIndex)
    words = shannon(firstIndex,index,p,words);
    data = words;
end

if(endIndex > index)
    words = shannon(index+1,endIndex,p,words);
    data = words;
end

else 
    data = words;
end

end

function inverted = varon(vector)
inverted = zeros(1,length(vector));
for i = 1:length(vector)
    inverted(i) = vector(length(vector) + 1 - i );
end
end


function codedImage = code(image,map)
codedImage = string(zeros(size(image)));
[Isize,Jsize] = size(image);
for i = 1:Isize
    for j = 1:Jsize
    image(i,j);
    codedImage(i,j) = map(image(i,j));
    end
end
end

function decodedImage = decode(codedImage,map)
decodedImage = zeros(size(codedImage));
[Isize,Jsize] = size(decodedImage);
for i = 1:Isize
    for j = 1:Jsize
    decodedImage(i,j) = map(codedImage(i,j));
    end
end
decodedImage = uint8(decodedImage)
end

function modulated = minModule(m,fc,Ts,fs)

% detecting Amc
if (m == 1 || m == 4)
    Amc = 1;
else
    Amc = -1;
end

% detecting Ams   
if(m == 1 || m == 2)
    Ams = 1;
else
    Ams = -1;
end
% sampling rate?!
t = 1/fs:1/fs:Ts;
modulated = Amc*sqrt(2/Ts).*cos(2*pi*fc*t) + Ams*sqrt(2/Ts).*sin(2*pi*fc*t);
end


% function repositoryMap = getRep(m,fc,Ts,fs)
% key = 1:1:4;
% valueSet = zeros(size(key));
% for i = 1:length(key)
%     valueSet(i) = module(m,fc,Ts,fs);
% end
% repositoryMap = containers.Map(key,valuseSet);
% end


%% todo 
function out = minDemodulateHilbert(signal,fc,Ts,fs)
signalHilbert = imag(hilbert(signal));
domain = 1/fs:1/fs:Ts;
AcMat = signal.*cos(2*pi*fc*domain) + signalHilbert.*sin(2*pi*fc*domain);
Ac = sign(AcMat(1,1));

AsMat = signal.*sin(2*pi*fc*domain) - signalHilbert.*cos(2*pi*fc*domain);
As = sign(AsMat(1,1));
out = [Ac,As];
end


function [out1,out2] = minDemodulateCorr(signal,fc,Ts,fs,phi)
sample_times = 1/fs:1/fs:Ts;
a = sqrt(2 / Ts);
n = ceil(length(signal) / 2);

s11 = a.*(cos(2*pi*fc*sample_times + phi)+sin(2*pi*fc*sample_times + phi));
s12 = a.*(cos(2*pi*fc*sample_times + phi)-sin(2*pi*fc*sample_times + phi));
s21 = a.*(-cos(2*pi*fc*sample_times + phi)+sin(2*pi*fc*sample_times + phi));
s22 = a.*(-cos(2*pi*fc*sample_times + phi)-sin(2*pi*fc*sample_times + phi));

c11 = sum(s11 .* signal);
c12 = sum(s12 .* signal);
c21 = sum(s21 .* signal);
c22 = sum(s22 .* signal);

c = [c11,c12,c21,c22];
[~, i] = max(c);

a_ap = xcorr(cos(2*pi*fc*sample_times), signal) ./ fs;
b_ap = xcorr(sin(2*pi*fc*sample_times), signal) ./ fs;
a_ap = a_ap(1, n) * 2 * a;
b_ap = b_ap(1, n) * 2 * a;

if i == 1
    Ac = 1; As = 1;
end
if i == 2
    Ac = 1; As = -1;
end
if i == 3
    Ac = -1; As = 1;
end
if i == 4
    Ac = -1; As = -1;
end
out1 = [Ac,As];
out2 = [a_ap,b_ap]; % -> Ac As aprox.
end

function m = detector(constellation)
Ac = constellation(1,1);
As = constellation(1,2);
m = 0;
if(Ac == 1 && As == 1)
    m = 1;
elseif (Ac == 1 && As == -1)
    m = 4;
elseif (Ac == -1 && As == 1)
    m = 2;
elseif(Ac == -1 && As == -1)
    m = 3;
end

end

function output = createStringBinary(mat,noise)
if(noise == 0)
    
output = "";
[Iindex,Jindex] = size(mat);
for i = 1:Iindex
    for j = 1:Jindex
        output = output + mat(i,j);
    end
end

else 
    output = "#";   
    [Iindex,Jindex] = size(mat);
    for i = 1:Iindex
        for j = 1:Jindex
            output = output + mat(i,j) + "#";
        end
    end
    
end
end

function output = extract(binaryString,map,Iindex,Jindex,noise)
if(noise == 0)
binaryChars = char(binaryString);
len = length(char(binaryString));
codedImage = zeros(Iindex,Jindex);
endindex = 1;
startindex = 1;
i = 1;
j = 1;
do = 1;
while do == 1
    substring = binaryChars(startindex:endindex);
    exist = isKey(map,substring);
    
    if(exist)
        codedImage(i,j) = map(substring);
        j = j+1;
        startindex = endindex + 1;
        endindex = startindex;
        if(endindex >len)
            break;
        end
    else
        endindex = endindex + 1;
        if(endindex > len)
            break;
        end
    end
    
    if(j > Jindex)
       i = i+1;
       j = 1;
       if(i > Iindex)
           break;
       end
    end
   
end
output = uint8(codedImage);

else
binaryChars = char(binaryString);
codedImage = zeros(Iindex,Jindex);
hashtags = find(binaryChars == '#');
jindex = 1;
iindex = 1;
for i = 1:length(hashtags)-1
    myString = string(binaryChars(hashtags(1,i)+1:hashtags(1,i+1)-1));
    if(isKey(map,myString))
        codedImage(iindex,jindex) = map(myString);
        jindex = jindex + 1;
    else
        codedImage(iindex,jindex) = nearest(myString,map);
        jindex = jindex + 1;
    end
    if(jindex > Jindex)
        jindex = 1;
        iindex = iindex + 1;
        if(iindex > Iindex)
            break;
        end
    end
end
output = uint8(codedImage);
end

end


function output = nearest(code,map)
keySet = keys(map);
len = length(keySet);
dists = zeros(1,len);
for i = 1:len
    dists(1,i) = distNotsynced(char(code),char(keySet(i)));
end
minIndex = find(dists == min(dists));
output = map(string(keySet(minIndex(1))));
end


function out = distNotsynced(code1,code2)
char1s = char(code1);
char2s = char(code2);
len1 = length(char1s);
len2 = length(char2s);

if(len1 >= len2)
    A = char(zeros(1,len1));
    A(1:len2) = char2s;
    A(len2+1:len1) = '0';
    B = char1s;
else
    A = char(zeros(1,len2));
    A(1:len1) = char1s;
    A(len1+1:len2) = '0';
    B = char2s;
end
out = distsynced(A,B);
end


function out = distsynced(A,B)
counter = 0;
for i = 1:length(A)
    if(A(i) ~= B(i) )
        counter = counter + 1;
    end
end
out = counter;
end

function out = syncer(chars)
if(chars(1) ~= '#')
    out = chars;
else
    chars(find(chars == '#')) = '';
    out = chars;
end
end

function  out = modulate(binaryString,fs,fc,Ts)
eachAnalSignalSize = Ts*fs;
startCursor = 1;
finalCursor = eachAnalSignalSize;
chars = char(binaryString);
chars = syncer(chars);
binaryString = string(chars);

len = length(chars);
if(mod(len,2) ~= 0)
    chars = char(binaryString + "0");
end
numberofSignals = (length(chars)/2);
sizeoftotalAnalSignal = numberofSignals * eachAnalSignalSize;
analSig = zeros(1,sizeoftotalAnalSignal);
for i = 0:numberofSignals-1        
    m = mFinder(chars(2*i+1:2*(i+1)));
    analSig(startCursor:finalCursor) = minModule(m,fc,Ts,fs);
    startCursor = startCursor + eachAnalSignalSize;
    finalCursor = finalCursor + eachAnalSignalSize;        
end
out = analSig;
end


function out = mFinder(chars)

% consider binary code as char1char2
if(chars(1) == '0' && chars(2) == '0')
    out = 1;
elseif(chars(1) == '0' && chars(2) == '1')
    out = 2;
elseif(chars(1) == '1' && chars(2) == '0')
    out = 3;
else
    out = 4;
end

end

function out = m2binary(m)
if( m == 1 )
    out = "00";
elseif( m == 2 )
    out = "01";
elseif( m == 3 )
    out = "10";
else % m == 4
    out = '11';
end
end

function out = demodule(signal,noisy,fs,fc,Ts,flags,phase,draw)

eachAnalSignalSize = Ts*fs;
startCursor = 1;
finalCursor = eachAnalSignalSize;
numberOf2bitSignal = length(signal)/eachAnalSignalSize;
rxSig = zeros(1,numberOf2bitSignal);
binaryCoded = char(zeros(1,2*numberOf2bitSignal));
    for i = 0:numberOf2bitSignal-1
        analSig = signal(startCursor:finalCursor);
        [out1,out2] = minDemodulateCorr(analSig,fc,Ts,fs,phase);
        rxSig(1,i+1) = out2(1,1) + 1j*out2(1,2);       
        binaryCoded(2*i+1:2*(i+1)) = char(m2binary(detector(out1)));
        startCursor = startCursor + eachAnalSignalSize;
        finalCursor = finalCursor + eachAnalSignalSize;       
    end
    
if(draw == 1) 
realLine = -2:0.01:2;
imagLine = 1j*(-2:0.01:2);

c = [1+1j 1-1j -1+1j -1-1j]
h = scatterplot(rxSig);
hold on
scatterplot(c,[],[],'r*',h)
hold on 
scatterplot(realLine,[],[],'g*',h)
scatterplot(imagLine,[],[],'g*',h)

grid

hold off
end

if(noisy == 1)
% todo add flags ...
output = "#"';
startIndex = 1;
for i = 1:length(flags)-1
    len = flags(i+1)-flags(i) - 2 ;
    endIndex = startIndex + len ;
    output = output + string(binaryCoded(startIndex:endIndex)) + "#";
    startIndex = endIndex + 1;
end
out = output;
return;
end
out = string(binaryCoded);
end

function out = Channel(signal,fc,fs,w)
fpassDown = fc - w/2;
fpassUp = fc + w/2;
out = bandpass(signal,[fpassDown fpassUp],fs);
end

function [out,Tsignal] = sdf(signal,fs)
[Tsignal,~] = standardFourierTransform(signal,fs);
out = Tsignal.*conj(Tsignal);
end

function [out,standardDomain] = standardFourierTransform(signal,fs)
n = length(signal);
fourier=fft(signal,n);
fourier = fourier/fs;
out = fftshift(fourier);
domain = linspace(-pi,pi,length(out));
standardDomain = domain*fs/(2*pi);
end

function boundedEnergy = calEnergy(powerSignal,center,range)
boundedEnergy = sum(powerSignal(center-range:center+range));
end



function bw = BW(signal,fs)
[~,finalDomain] = standardFourierTransform(signal,fs);
[audioDataSdf,~] = sdf(signal,fs);
tempCenter = find(finalDomain >= 10000);
centerFreq = tempCenter(1);
zeroTemp = find(finalDomain >= 0);
zero = zeroTemp(1);
powerAbs = abs(audioDataSdf);
totalEnergy = sum(powerAbs(zero:end));
power99 = totalEnergy*0.99;
freqIterator = 100;% 1kHz
boundedEnergy = calEnergy(powerAbs,centerFreq,freqIterator);
while boundedEnergy < power99
    freqIterator = freqIterator + 1000;
    boundedEnergy = calEnergy(powerAbs,centerFreq,freqIterator);
end
startBound = centerFreq - freqIterator;
endBound = centerFreq + freqIterator;
normalizedStartBound = finalDomain(startBound:endBound);
bw = (vpa(normalizedStartBound(end)) - vpa(normalizedStartBound(1))) % it is bandPass there is no need to /2;
end

