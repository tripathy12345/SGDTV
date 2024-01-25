function psd_fi=MAE_PSD_all_bands(x1, receeg, Fs)
Nfft=Fs; 
window=kaiser(64,3);  
noverlap=32;   
% [Pxx,fp]=pwelch(x2,window,noverlap,Nfft,Fs);
[Pf,fp]=pwelch(receeg,window,noverlap,Nfft,Fs);
[Pref,fp]=pwelch(x1,window,noverlap,Nfft,Fs);

% Pn=10*log10(Pxx);
Pf=10*log10(Pf);
Pref=10*log10(Pref);

% figure    
% plot(fp,10*log10(Pn),fp,10*log10(Pf),'-r', fp,10*log10(Pref),'-k');  
% xlim([0 30])
% xlabel('Frequency (Hz)');ylabel('Power Spectrum (dB)');
% title('Power Spectrum of delta Rhythm'); grid on
% legend('Noisy signal','MSGTV filter output', 'Reference');

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
res=[MAE_D MAE_T MAE_A MAE_B];
% psd_fi=sum(res);
psd_fi=MAE_A+ MAE_B;
end