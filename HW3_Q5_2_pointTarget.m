%% two parallel receive
pointTarget_data=pointTarget.data(80:end,:,:);
time=[0:1:size(pointTarget_data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));
timeArray=[0:1/(pointTarget.samplingRateMHz*(10^6)):(size(pointTarget_data,1)-1)/(pointTarget.samplingRateMHz*(10^6))]';
timeArray2=repmat(timeArray,[1,128]);

for zz=1:length(pointTarget_data)
    zf(zz,1)=(time(zz)*1540)/2;
end 

%%
%pos pitch 
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_pos_1(yy,bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        xf_pos_1(yy,bb+64.5)=(.5)*((pointTarget.elementSpacingMM)/1000);;
        diag_dist_pos_1(yy,bb+64.5)=sqrt(zf(yy)^2 + ((xe_pos_1(yy,bb+64.5)-xf_pos_1(yy,bb+64.5))^2));
        time_diag_pos_1(yy,bb+64.5)=diag_dist_pos_1(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_pos_1(yy,bb+64.5)=time_diag_pos_1(yy,bb+64.5)-time_diag_pos_1(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_pos_1(dd,:)=timeArray2(dd,:)+time_delay_pos_1(dd,:);
end

for hh=1:128
    temp_pos_1=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:2:128)),time_withDelays_pos_1(:,hh),'linear',0);
    reshaped_interp_pos_1(:,hh,:)=reshape(temp_pos_1,[2353,1,64]);
end 


%%
%neg pitch 
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_neg_1(yy,bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        xf_neg_1(yy,bb+64.5)=(-.5)*((pointTarget.elementSpacingMM)/1000);
        diag_dist_neg_1(yy,bb+64.5)=sqrt(zf(yy)^2 + ((xe_neg_1(yy,bb+64.5)-xf_neg_1(yy,bb+64.5))^2));
        time_diag_neg_1(yy,bb+64.5)=diag_dist_neg_1(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_neg_1(yy,bb+64.5)=time_diag_neg_1(yy,bb+64.5)-time_diag_neg_1(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_neg_1(dd,:)=timeArray2(dd,:)+time_delay_neg_1(dd,:);
end

for hh=1:128
    temp_neg_1=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:2:128)),time_withDelays_neg_1(:,hh),'linear',0);
    reshaped_interp_neg_1(:,hh,:)=reshape(temp_neg_1,[2353,1,64]);
end 

%% 
interp_2=zeros(2353,128,128);

for mm=1:64
   interp_2(:,:,2*mm-1)=reshaped_interp_neg_1(:,:,mm);
   interp_2(:,:,2*mm)=reshaped_interp_pos_1(:,:,mm);
end 

 for jj=1:128
     for kk=1:2353
            summed_interp_2(kk,jj)=sum(interp_2(kk,:,jj));
     end
end

figure;
imagesc(20*log10(abs(hilbert(summed_interp_2))),[30,80])
title('2 parallel receive beams')
colormap('gray')
