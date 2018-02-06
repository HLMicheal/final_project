function[]=findfb()
T=0.6;  % sampling period
Fs=1e6; % sampling freq

    % generate FMCW
f1=2.4e9; % min freq
f2=2.5e9; % max freq
B=f2-f1;  % bandwidth
% N=T*Fs+1;
N=T*Fs;   % number of steps
D=1.5e-8;   % changable, time delay
for i=1:1:N
t=(i-1)/Fs;
A=f(t);
y=cos(2*pi*A*t);
VTx=y/2;
    % mix the transmitted and received signals
C=f(t-D);
VRx=0.5*cos(2*pi*C*(t-D));
if abs(VTx-VRx)>0.001
    output=1;
end
mixed(i)=VTx*VRx;   % mixed is the received signal, in other words, the only data I have
end
range=3e8*D/2;  % calculate the distance of the object

    % select a time to calculate the beating frequency
S=fft(mixed,N);
ST=fftshift(S); % Fourier transform
% fshift=(-(N-1)/2:(N-1)/2)*(Fs/N);
fshift=(-N/2:N/2)*(Fs/N);  % x dimension 
ff=abs(fshift)>16000;   
ST(ff)=0;    % limit of the max range that can be detected

    % calculate the range
[Y1,I1]=max(abs(ST));  % find the max value of ST
fb=abs(fshift(I1)); % the corresponding freq
if fb>=36
    ft=abs(fshift)>fb/2;
    ST(ft)=0;
    [Y2, I2]=max(abs(ST));
    fb=abs(fshift(I2));
end
requiredfb=range*2*B/(3e8*0.04);
R=3e8*0.04*fb./(2*B);

figure(2)
plot(fshift,abs(ST))
xlim([-fb-100 fb+100])
grid on;
fprintf(' fb:        %g \n required fb: %g\n actual range:   %g\n measured range: %g\n',fb,requiredfb,range,R);
end

function [ y ] = f( t )
while (t<0)
    t=t+0.04;
end
while (t>0.04)
    t=t-0.04;
end
f1=2.4e9;
f2=2.5e9;
B=f2-f1;
alpha=B/0.04;
y=f1+alpha.*t;
end
