%% two parallel receive
anecoicCyst_data=anecoicCyst.data(80:end,:,:);
time=[0:1:size(anecoicCyst_data,1)-1]*(1/((anecoicCyst.samplingRateMHz)*(10^6)));
timeArray=[0:1/(anecoicCyst.samplingRateMHz*(10^6)):(size(anecoicCyst_data,1)-1)/(anecoicCyst.samplingRateMHz*(10^6))]';
timeArray2=repmat(timeArray,[1,128]);

for zz=1:length(anecoicCyst_data)
    zf(zz,1)=(time(zz)*1540)/2;
end 

%%
%pos pitch 
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_cont_pos(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos(yy,bb+64.5)=xe_cont_pos(bb+64.5)*(.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos(yy,bb+64.5))^2)-((xf_pos(yy,bb+64.5))^2)));
        time_diag_cont_pos(yy,bb+64.5)=diag_dist_cont_pos(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos(yy,bb+64.5)=time_diag_cont_pos(yy,bb+64.5)-time_diag_cont_pos(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos(dd,:)=timeArray2(dd,:)+time_delay_cont_pos(dd,:);
end

for hh=1:128
    temp_cont_pos=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:2:128)),time_withDelays_cont_pos(:,hh),'linear',0);
    reshaped_interp_cont_pos(:,hh,:)=reshape(temp_cont_pos,[2353,1,64]);
end 

for jj=1:64
    for kk=1:2353
            zone_interp_cont_2_pos(kk,jj)=sum(reshaped_interp_cont_pos(kk,:,jj));
    end
end

%%
%neg pitch
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_cont_neg(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg(yy,bb+64.5)=xe_cont_neg(bb+64.5)*(-.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg(yy,bb+64.5))^2)-((xf_neg(yy,bb+64.5))^2)));
        time_diag_cont_neg(yy,bb+64.5)=diag_dist_cont_neg(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg(yy,bb+64.5)=time_diag_cont_neg(yy,bb+64.5)-time_diag_cont_neg(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg(dd,:)=timeArray2(dd,:)+time_delay_cont_neg(dd,:);
end

for hh=1:128
    temp_cont_neg=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:2:128)),time_withDelays_cont_neg(:,hh),'linear',0);
    reshaped_interp_cont_neg(:,hh,:)=reshape(temp_cont_neg,[2353,1,64]);
end 

for jj=1:64
    for kk=1:2353
            zone_interp_cont_2_neg(kk,jj)=sum(reshaped_interp_cont_neg(kk,:,jj));
    end
end

%% 
zone_interp_cont_2=zeros(2353,128);

for mm=1:64
    zone_interp_cont_2(:,2*mm-1)=zone_interp_cont_2_pos(:,mm);
    zone_interp_cont_2(:,2*mm)=zone_interp_cont_2_neg(:,mm);
end 

figure;
imagesc(20*log10(abs(hilbert(zone_interp_cont_2))),[30,80])
title('2 parallel receive beams')
colormap('gray')
