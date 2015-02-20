for zf=[.015 .03 .045 .06 .075 ]

time=[0:1:size(pointTarget.data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));

rate_upsample=((pointTarget.samplingRateMHz)*(10^6))*4;

time_upsample=[0:1/4:size(pointTarget.data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));

figure;
for aa=1:128
    interpolation(:,:,aa)=interp1(time,pointTarget.data(:,:,aa),time_upsample,'linear');
    imagesc(interpolation(:,:,aa))
    colormap('gray')
end

for bb=0:63
    xe_1(bb+1)=((pointTarget.elementSpacingMM)/1000)*(64-bb);
    diag_dist_1(bb+1)=sqrt(zf^2 + (xe_1(bb+1))^2);
    time_diag_1(bb+1)=diag_dist_1(bb+1)/1540; 
end

for cc=0:63
    time_delay_1(cc+1)=time_diag_1(cc+1)-time_diag_1(64);
end 

for dd=1:64
    xe_2(dd)=((pointTarget.elementSpacingMM)/1000)*dd;
    diag_dist_2(dd)=sqrt(zf^2 + (xe_2(dd))^2);
    time_diag_2(dd)=diag_dist_2(dd)/1540;
    time_delay_2(dd)=time_diag_2(dd)-time_diag_2(1);
end

xe=[xe_1 xe_2];
diag_dist=[diag_dist_1 diag_dist_2];
time_diag=[time_diag_1 time_diag_2];
time_delay=[time_delay_1 time_delay_2]; 

depth_range_1=[0:zf/63:zf]
depth_range_2=[zf:-zf/63:0]
depth_range=[depth_range_1 depth_range_2]

for ee=1:length(depth_range)
   timeArray(ee)=depth_range(ee)/rate_upsample
end

for ff=1:length(depth_range)
timeArray_withDelay(ff)=timeArray(ff)+time_delay(ff)
end 

for gg=1:128
    delayedChannel=interp1(timeArray_withDelay,pointTarget.data(:,:,gg),timeArray,'linear')
    imagesc(interpolation(:,:,gg))
    colormap('gray')
end

for 1:128
    time=time




end