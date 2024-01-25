clc;
clear all;
close all;
tic
load reference.mat;
load contaminated.mat;
%%%x1 is noise free and x2 is corrupted signal.
Fs=200;
f=x2';
N=length(f);
params.SamplingRate=Fs;

L=6;
frequency=(Fs./(2.^(L:-1:1)));
 % frequency=[4 8 13 30 75]; %%%%For rhythm specific case
subresf=1;
boundaries=(frequency*(2*pi))/Fs;

w1=frequency(1:end-1);
w2=frequency(2:end);
gamma=(w2-w1)./(w2+w1);


ff=fft(f);
% We build the corresponding filter bank
mfb=EWT_Meyer_FilterBank(boundaries,length(ff));

% We filter the signal to extract each subband
ewt=cell(length(mfb),1);
for k=1:length(mfb)
    mm=real(ifft(conj(mfb{k}).*ff));
    ewt{k}=mm;
    modes(k,:)=mm;
end

Bound=1;
xxx=(linspace(0,1,round(length(mfb{1,1}))))*Fs;
for i=1:size(mfb)
plot(xxx,mfb{i,1})
hold on
end
xlim([0 round(Fs/2)])
ylim([0 2])
title('DBPEWT filter bank')
div=1;
Show_EWT_Boundaries(abs(fft(f)),boundaries,div,params.SamplingRate);

% figure
% % subplot(521)
% % plot(f);
% 
% subplot(611)
% plot(ewt{1,1});
% 
% subplot(612)
% plot(ewt{2,1});
% subplot(613)
% plot(ewt{3,1});
% subplot(614)
% plot(ewt{4,1});
% subplot(615)
% plot(ewt{5,1});

scale1=ewt{1,1};
scale2=ewt{2,1};
scale3=ewt{3,1};
scale4=ewt{4,1};
scale5=ewt{5,1};
scale6=ewt{6,1};
order=2;
alpha=0.15;
window=round(length(f)/(alpha*Fs)); 
if (rem(window,2)==0)
    window=window+1;
end
xsg = sgolayfilt(scale1,order,window);
lam = 0.15;                         % lam: regularization parameter
Nit = 50;                        % Nit: number of iterations
[xaf, cost] = tvd_mm(xsg, lam, Nit);

deltaf=scale1-xaf;
figure
plot(scale1)
hold on
plot(deltaf)

receeg=deltaf+scale2+scale3+scale4+scale5+scale6;
MIbf=mi(x2',receeg) %%%%filtered and contaminated
MIaf=mi(x1', receeg) %%%filtered and reference
R_bf = corrcoef(x2', receeg); %%%contaminated and filtered
rowbf=R_bf(1,2);
R_af = corrcoef(x1', receeg); %%%%reference and filtered
rowaf=R_af(1,2);

SARbf= 10*log(std(x2)./std(x2'-receeg));
SARaf= 10*log(std(x1)./std(x1'-receeg));


toc
Nfft=Fs; 
window=kaiser(100,3);  
noverlap=50;   
[Pxx,fp]=pwelch(x2,window,noverlap,Nfft,Fs);
[Pf,fp]=pwelch(receeg,window,noverlap,Nfft,Fs);
[Pref,fp]=pwelch(x1,window,noverlap,Nfft,Fs);

Pn=10*log10(Pxx);
Pf=10*log10(Pf);
Pref=10*log10(Pref);

figure    
plot(fp,10*log10(Pn),fp,10*log10(Pf),'-r', fp,10*log10(Pref),'-k');  
xlim([0 30])
xlabel('Frequency (Hz)');ylabel('Power Spectrum (dB)');
title('Power Spectrum of delta Rhythm'); grid on
legend('Noisy signal','MSGTV filter output', 'Reference');

for i=1:5
MAE_d(i)=abs(Pref(i)-Pf(i))./(5-0);
end
MAE_D=sum(MAE_d);

for j=6:9
MAE_th(i)=abs(Pref(i)-Pf(i))./(4-0);
end
MAE_T=sum(MAE_th);

for k=10:15
MAE_a(k)=abs(Pref(k)-Pf(k))./(6-0);
end
MAE_A=sum(MAE_a);
for kk=16:31
MAE_be(kk)=abs(Pref(kk)-Pf(kk))./(18-0);
end
MAE_B=sum(MAE_be);
res=[MAE_D MAE_T MAE_A MAE_B]

