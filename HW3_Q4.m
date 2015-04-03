figure;
plot(20*log10(abs(fftshift(fft(zone_interp_cont(:,65))))));
title('Mag of FT of RF beam 65')
 
for gg=1:128
    fourier_rf(:,gg)=20*log10(abs(fftshift(fft(zone_interp_cont(:,gg)))))';
end 

mean_fourier_rf=mean(fourier_rf,2);

figure; 
plot(mean_fourier_rf)
title('Avg Mag of FT of all RF Beams')


figure;
b=fir1(10,228/2353,'bandpass')
c=filtfilt(b,1,20*log10(abs(fftshift(fft(zone_interp_cont(:,65))))));
plot(c)
title('filtered RF data beam 65')

figure;
d=fir1(10,228/2353,'bandpass')
e=filtfilt(d,1,mean_fourier_rf);
plot(e)
title('filtered RF data all beams')

figure;
f=fir1(10,228/2353,'bandpass')
g=filtfilt(f,1,20*log10(abs(hilbert(zone_interp_cont))));
plot(g)
imagesc(lat_array,axial_array,g)
colormap('gray')
axis image
title('filtered cyst')
