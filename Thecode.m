%dlmread('pressure_experiment(1).dat');
%dlmread('pressure_raw_cfd(1).dat');
time1_t=pressurerawcfd1(:,1);    % readfing cfd data
P1_t=pressurerawcfd1(:,2);       %Reading cfd data

f2=pressureexperiment1(:,1);     %rading experimental data
P2=pressureexperiment1(:,2);     %reading experimental data

Num=length(P1_t);            % number of data points in given data
v=mod(Num,2);                % remainder after dividing number by 2, either 1 or 0 
if v==1                      % basically if number not divisible by 2, substract 1 from it
Num=Num-1;                   % basically if number not divisible by 2, substract 1 from it
end                          % basically if number not divisible by 2, substract 1 from it
% Pr1_w1=[1:Num];              % making new array of 12345678.....num      
% time1=[1:Num];               % making new array of 12345678.....num
  for i=1:Num
   Pr1_w1(i)=P1_t(i) ;     %copying the array of p1_t to pr1_w1
   time1(i)=time1_t(i);       %copying the array from time1_t to time1 
  end
  Fs=Num/(time1(Num)-time1(1)); % frequency = total number of data points / time last - time first, total length of data by total time f=1/time
  Pr1_w2=Pr1_w1;              % copying the data from one array to other
  Pm=mean(Pr1_w1);            % mean value, one value of mean of pr1_w1
 Pr1_w1=Pr1_w1-Pm;           % deviation from mean, not absolute value
 Pr1_w2=Pr1_w2-Pm;           % deviation from the mean, the same array as the pr_1
%  Windows
 Y_w1=transpose(hamming(Num)); %hamming Window  % transpose- interchange row and columns and hamming of numbers, which is, the start number and end num is 0.08-1-0.08, and the rest are equally spaced between
 Pr1_w1= Pr1_w1.*Y_w1     % hamming Window     %hamming window value is multiplied value is multiplied to the deviation from mean
 Y_w2=transpose(hann(Num)); %hanning Window    % transpose is done to make it a row vector from column vector, hann is window (0-1-0) 
 Pr1_w2= Pr1_w2.*Y_w2; %hann Window % the hanning window value is multiplied to the deviation from mean       
 Pr1_w1=Pr1_w1+Pm;  % mean is added to the windowed value
 Pr1_w2=Pr1_w2+Pm  % mean is added to the windowed value
% FFT
 Y_w1=fft(Pr1_w1);           %fft applied to hamming
 Y_w2=fft(Pr1_w2);           %fft applied to hann
 Pw2_w1=abs(Y_w1/Num);       % absolute value of fft divided by the total number
 Pw2_w2=abs(Y_w2/Num);       % absolute value of fft divided by the total number
 Pw1_w1=Pw2_w1(1:Num/2+1);   % from 1 to first half + 1 values
 Pw1_w2=Pw2_w2(1:Num/2+1);   % from 1 to first half + 1 values
 Pw1_w1(2:end-1)=2*Pw1_w1(2:end-1);   %from 2 to end - 1 
 Pw1_w2(2:end-1)=2*Pw1_w2(2:end-1);   %from 2 to end - 1
% Pw1_w1=Pw1_w1.*2;               %hamming Window Power Correction
%Pw1_w2=Pw1_w2.*2;                %hanning window power corection
 f=Fs*(0:(Num/2))/Num;              %freq multiplied by  0-0.5
 SPL_w1=10*log10(Pw1_w1.^2/(4*(10^(-10)))); %converting to SPL hamming
 SPL_w2=10*log10(Pw1_w2.^2/(4*(10^(-10)))); %converting to SPL hann
 SPL=(SPL_w1+SPL_w2); %Ensambling Average of both the windows
 SPL=SPL./2; %ensambling Average
 plot(f,SPL) % frequency vs SPL
%xlim([1500 4500])
hold on
%plot(f2,P2)
hold off
f=transpose(f);
SPL=transpose(SPL);
SS=[f,SPL];
dlmwrite('avg_spl.dat',SS);
