%Reading data 
PCFD= pressurerawcfd1(:,2); % CFD pressure data 
TCFD= pressurerawcfd1(:,1); % CFD time data 
fexp=pressureexperiment1(:,1); %rading experimental data 
SPLexp=pressureexperiment1(:,2); %reading experimental data 
Fpass=3000 ; %according to the paper 
newFsig=6000 ; %nyquest frequency 
dt=pressurerawcfd1(2,1)-pressurerawcfd1(1,1); 
Fsig=1/dt; 
%low pass filter 
P3 = lowpass(PCFD,Fpass,Fsig) ; %low pass filter 
%resampling 
[p,q]=rat(newFsig/Fsig); 
P22=resample(P3,p,q); 
t1=(0:(length(P22)-1))*q/(p*(1/dt)); 
figure(100); 
plot(TCFD,P3,t1,P22); 
legend("Original","resample") 
title("resample") 
num=length(PCFD); 
num=num-1; 
for i=1:num 
Phamm(i)=PCFD(i); 
Pnowin(i)=PCFD(i); 
TCFD(i)=TCFD(i); 
end 
Fs=num/TCFD(num)-TCFD(1); %frequency 
Pm=mean(Phamm); % mean pressure 
Pdev=Phamm-Pm; %deviation of pressure from mean 
% %Pdev=Phann-Pm; %deviation of pressure from mean 
%number array 
n=num; 
array=[1:n] ; % create your array 
serialno=transpose(array); 
%Windowing 
tuk=transpose(tukeywin(num)); %Tukeywin 
ham=transpose(hamming(num)); %hamm window 0.08-1-0.08 
han=transpose(hann(num)); %hann window 0-1-0 
blkman=transpose(blackman(num)); %blackman window 
cheby=transpose(chebwin(num)); %chebyshev window 
bart=transpose(bartlett(num)); %bartlett window 0-1-0 
Pdevtuk=tuk.*Pdev ; % tukey window 
Pdevhamm=ham.*Pdev ; % hamming windowing on the pressure data 
Pdevhann=han.*Pdev; % hanning window on the pressure data 
Pdevblkman=blkman.*Pdev; % blkman window on the pressure data 
Pdevcheby=cheby.*Pdev; %cheby window 
Pdevbart=bart.*Pdev; %bartlett window 
Pw1=Pm+Pdevtuk; % pressure mean + the tuk window 
Pw2=Pm+Pdevhamm; % pressure mean + hamm window 
Pw3=Pm+Pdevhann; % pressure mean + hann window 
Pw4=Pm+Pdevblkman; % blackman 
Pw5=Pm+Pdevcheby; %cheby 
Pw6=Pm+Pdevbart; %bart 
% 
%FFT 
fftw0=fft(Pnowin); 
fftw1=fft(Pw1); %fft applied to tuk 
fftw2=fft(Pw2); %fft applied to hamm 
fftw3=fft(Pw3); 
fftw4=fft(Pw4); 
fftw5=fft(Pw5); 
fftw6=fft(Pw6); 
Pfftw0=abs(fftw0/num); 
Pfftw1=abs(fftw1/num); % absolute value of fft divided by the total number 
Pfftw2=abs(fftw2/num); % absolute value of fft divided by the total number 
Pfftw3=abs(fftw3/num); 
Pfftw4=abs(fftw4/num); 
Pfftw5=abs(fftw5/num); 
Pfftw6=abs(fftw6/num); 
P2fftw0=Pfftw0(1:num/2); 
P2fftw1=Pfftw1(1:num/2); %from 1 to first half 
P2fftw2=Pfftw2(1:num/2); %from 1 to first half 
P2fftw3=Pfftw3(1:num/2); 
P2fftw4=Pfftw4(1:num/2); 
P2fftw5=Pfftw5(1:num/2); 
P2fftw6=Pfftw6(1:num/2); 
f=Fs*(0:(num/2-1))/num; %freq multiplied by 0-0.5 to reduce leakage 
Pref=0.00002; % preference which is given 
SPLnowin=20*log10(P2fftw0/0.00002); 
SPLw1=20*log10(P2fftw1/0.00002); %converting to SPL 
SPLw2=20*log10(P2fftw2/0.00002); %converting to SPL 
SPLw3=20*log10(P2fftw3/0.00002); 
SPLw4=20*log10(P2fftw4/0.00002); 
SPLw5=20*log10(P2fftw5/0.00002); 
SPLw6=20*log10(P2fftw6/0.00002); 
SPLw1 = sgolayfilt(SPLw1,5,17); %sgolay filter with 5 and 17 averaging 
SPLw2 = sgolayfilt(SPLw2,5,17); 
SPLw3 = sgolayfilt(SPLw3,5,17); 
SPLw4 = sgolayfilt(SPLw4,5,17); 
SPLw5 = sgolayfilt(SPLw5,5,17); 
SPLw6 = sgolayfilt(SPLw6,5,17); 
figure(10) 
plot(f,SPLw1,'b') 
xlim([1500 4500]) 
title("All windows sgolay filter applied") 
xlabel("SPL (dB)") 
ylabel("frequency") 
hold on 
plot(f,SPLw2,'m') 
xlim([1500 4500]) 
hold on 
plot(f,SPLw3,'r') 
xlim([1500 4500]) 
hold on 
plot(f,SPLw4,'g') 
xlim([1500 4500]) 
hold on 
plot(f,SPLw5,'y') 
xlim([1500 4500]) 
hold on 
plot(f,SPLw6,'c') 
xlim([1500 4500]) 
hold on 
plot(fexp,SPLexp,'--k') 
legend("Tukey window","Hamming window", "Hann window","Black man window","Chebyshev window","Bartlett window","Experimental results") 
figure(2) 
plot(f,SPLw2,'b') % frequency vs SPL 
title("CFD data hamming window") 
xlabel("SPL (dB)") 
ylabel("frequency") 
legend("CFD data") 

%Plotting 
%Cd plotting 
figure (1) 
x10=coarse2cd111(:,1); 
y10=coarse2cd111(:,2); 
title("Second order scheme Cd calculation") 
xlabel=("Iterations") 
ylabel=("Cd") 
plot(x10,y10) 
hold on 
x11=medium2cd111(:,1); 
y11=medium2cd111(:,2); 
plot(x11,y11) 
hold on 
x12=fine2cd111(:,1); 
y12=fine2cd111(:,2); 
plot(x12,y12) 
hold on 
x13=finer2cd111(:,1); 
y13=finer2cd111(:,2); 
plot(x13,y13) 
hold on 
legend("coarse","medium","fine","finer") 
%Cl plotting 
x20=coarse2cl111(:,1); 
y66=coarse2cl111(:,2); 
title("Second order scheme Cl calculation") 
xlabel=("Iterations") 
ylabel=("Cd") 
plot(x20,y66) 
hold on 
x31=medium2cl111(:,1); 
y31=medium2cl111(:,2); 
plot(x31,y31) 
hold on 
x32=fine2cl111(:,1); 
y32=fine2cl111(:,2); 
plot(x32,y32) 
hold on 
x33=finer2cl111(:,1); 
y3=finer2cl111(:,2); 
plot(x33,y3) 
hold on 
legend("coarse","medium","fine","finer") 
%Cp plotting 
x66=coarse2cp111(:,1); 
y66=coarse2cp111(:,2); 
y67=-1*y66; 
title("Second order scheme Cp data") 
xlabel=("Iterations") 
ylabel=("Cd") 
scatter(x66,y67) 
hold on 
x31=medium2cp111(:,1); 
y31=medium2cp111(:,2); 
y31=-1*y31 
scatter(x31,y31) 
hold on 
x32=fine2cp111(:,1); 
y32=fine2cp111(:,2); 
y38=-1*y32 
scatter(x32,y38) 
hold on 
x33=finer2cp111(:,1); 
y33=finer2cp111(:,2); 
y33=-1*y33 
scatter(x33,y33) 
hold on 
x55=airfoilagarddata(:,1); 
y55=airfoilagarddata(:,2); 
scatter(x55,y55) 
hold on 
legend("coarse","medium","fine","finer","experimental") 
% %Cp plotting 
% x1=coarse2cp(:,1); 
% y1=coarse2cp(:,2); 
% plot (x1,y1) 
% hold on 
% x2=medium2cp(:,1); 
% y2=medium2cp(:,2); 
% plot (x2,y2) 
% hold on 
% x3=fine2cp(:,1); 
% y3=fine2cp(:,2); 
% plot (x3,y3) 
% hold on 
% x4=airfoilagarddata(:,1); 
% y4=airfoilagarddata(:,2); 
% plot (x4,y4) 
% hold on 
% legend("coarse","medium","fine","experimental") 