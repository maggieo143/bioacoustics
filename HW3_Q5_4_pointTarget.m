%% 4 parallel receive
pointTarget_data=pointTarget.data(80:end,:,:);
time=[0:1:size(pointTarget_data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));
timeArray=[0:1/(pointTarget.samplingRateMHz*(10^6)):(size(pointTarget_data,1)-1)/(pointTarget.samplingRateMHz*(10^6))]';
timeArray2=repmat(timeArray,[1,128]);

for zz=1:length(pointTarget_data)
    zf(zz,1)=(time(zz)*1540)/2;
end 
%%
% +1/2 pitch 
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
    temp_pos_1=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:4:128)),time_withDelays_pos_1(:,hh),'linear',0);
    reshaped_interp_pos_1(:,hh,:)=reshape(temp_pos_1,[2353,1,32]);
end 

%% 
% +3/2 pitch 
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_pos_2(yy,bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        xf_pos_2(yy,bb+64.5)=(1.5)*((pointTarget.elementSpacingMM)/1000);;
        diag_dist_pos_2(yy,bb+64.5)=sqrt(zf(yy)^2 + ((xe_pos_2(yy,bb+64.5)-xf_pos_2(yy,bb+64.5))^2));
        time_diag_pos_2(yy,bb+64.5)=diag_dist_pos_2(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_pos_2(yy,bb+64.5)=time_diag_pos_2(yy,bb+64.5)-time_diag_pos_2(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_pos_2(dd,:)=timeArray2(dd,:)+time_delay_pos_2(dd,:);
end

for hh=1:128
    temp_pos_2=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:4:128)),time_withDelays_pos_2(:,hh),'linear',0);
    reshaped_interp_pos_2(:,hh,:)=reshape(temp_pos_2,[2353,1,32]);
end 


%% four parallel receive
% -1/2 pitch 
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
    temp_neg_1=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:4:128)),time_withDelays_neg_1(:,hh),'linear',0);
    reshaped_interp_neg_1(:,hh,:)=reshape(temp_neg_1,[2353,1,32]);
end 

%% 
% -3/2 pitch 
for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_neg_2(yy,bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        xf_neg_2(yy,bb+64.5)=(-1.5)*((pointTarget.elementSpacingMM)/1000);
        diag_dist_neg_2(yy,bb+64.5)=sqrt(zf(yy)^2 + ((xe_neg_2(yy,bb+64.5)-xf_neg_2(yy,bb+64.5))^2));
        time_diag_neg_2(yy,bb+64.5)=diag_dist_neg_2(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_neg_2(yy,bb+64.5)=time_diag_neg_2(yy,bb+64.5)-time_diag_neg_2(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_neg_2(dd,:)=timeArray2(dd,:)+time_delay_neg_2(dd,:);
end

for hh=1:128
    temp_neg_2=interp1(timeArray2(:,hh),squeeze(pointTarget_data(:,hh,1:4:128)),time_withDelays_neg_2(:,hh),'linear',0);
    reshaped_interp_neg_2(:,hh,:)=reshape(temp_neg_2,[2353,1,32]);
end 
%%
interp_4=zeros(2353,128,128);

for mm=0:31
   interp_4(:,:,4*mm+1)=reshaped_interp_neg_2(:,:,mm+1);
   interp_4(:,:,4*mm+2)=reshaped_interp_neg_1(:,:,mm+1);
   interp_4(:,:,4*mm+3)=reshaped_interp_pos_1(:,:,mm+1);
   interp_4(:,:,4*mm+4)=reshaped_interp_pos_2(:,:,mm+1);
end 

 for jj=1:128
     for kk=1:2353
            summed_interp_4(kk,jj)=sum(interp_4(kk,:,jj));
     end
end

figure;
imagesc(20*log10(abs(hilbert(summed_interp_4))),[30,80])
title('4 parallel receive beams')
colormap('gray')
