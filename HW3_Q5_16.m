%% 16 parallel receive
anecoicCyst_data=anecoicCyst.data(80:end,:,:);
time=[0:1:size(anecoicCyst_data,1)-1]*(1/((anecoicCyst.samplingRateMHz)*(10^6)));
timeArray=[0:1/(anecoicCyst.samplingRateMHz*(10^6)):(size(anecoicCyst_data,1)-1)/(anecoicCyst.samplingRateMHz*(10^6))]';
timeArray2=repmat(timeArray,[1,128]);

for zz=1:length(anecoicCyst_data)
    zf(zz,1)=(time(zz)*1540)/2;
end
%%
%+15/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_1(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_1(yy,bb+64.5)=xe_cont_pos_1(bb+64.5)*(7.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_1(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_1(yy,bb+64.5))^2)-((xf_pos_1(yy,bb+64.5))^2)));
        time_diag_cont_pos_1(yy,bb+64.5)=diag_dist_cont_pos_1(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_1(yy,bb+64.5)=time_diag_cont_pos_1(yy,bb+64.5)-time_diag_cont_pos_1(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_1(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_1(dd,:);
end

for hh=1:128
    temp_cont_pos_1=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_1(:,hh),'linear',0);
    reshaped_interp_cont_pos_1(:,hh,:)=reshape(temp_cont_pos_1,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_1(kk,jj)=sum(reshaped_interp_cont_pos_1(kk,:,jj));
    end
end
%%
%+13/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_2(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_2(yy,bb+64.5)=xe_cont_pos_2(bb+64.5)*(6.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_2(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_2(yy,bb+64.5))^2)-((xf_pos_2(yy,bb+64.5))^2)));
        time_diag_cont_pos_2(yy,bb+64.5)=diag_dist_cont_pos_2(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_2(yy,bb+64.5)=time_diag_cont_pos_2(yy,bb+64.5)-time_diag_cont_pos_2(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_2(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_2(dd,:);
end

for hh=1:128
    temp_cont_pos_2=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_2(:,hh),'linear',0);
    reshaped_interp_cont_pos_2(:,hh,:)=reshape(temp_cont_pos_2,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_2(kk,jj)=sum(reshaped_interp_cont_pos_2(kk,:,jj));
    end
end
%%
%+11/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_3(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_3(yy,bb+64.5)=xe_cont_pos_3(bb+64.5)*(5.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_3(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_3(yy,bb+64.5))^2)-((xf_pos_3(yy,bb+64.5))^2)));
        time_diag_cont_pos_3(yy,bb+64.5)=diag_dist_cont_pos_3(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_3(yy,bb+64.5)=time_diag_cont_pos_3(yy,bb+64.5)-time_diag_cont_pos_3(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_3(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_3(dd,:);
end

for hh=1:128
    temp_cont_pos_3=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_3(:,hh),'linear',0);
    reshaped_interp_cont_pos_3(:,hh,:)=reshape(temp_cont_pos_3,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_3(kk,jj)=sum(reshaped_interp_cont_pos_3(kk,:,jj));
    end
end
%%
%+9/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_4(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_4(yy,bb+64.5)=xe_cont_pos_4(bb+64.5)*(4.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_4(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_4(yy,bb+64.5))^2)-((xf_pos_4(yy,bb+64.5))^2)));
        time_diag_cont_pos_4(yy,bb+64.5)=diag_dist_cont_pos_4(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_4(yy,bb+64.5)=time_diag_cont_pos_4(yy,bb+64.5)-time_diag_cont_pos_4(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_4(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_4(dd,:);
end

for hh=1:128
    temp_cont_pos_4=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_4(:,hh),'linear',0);
    reshaped_interp_cont_pos_4(:,hh,:)=reshape(temp_cont_pos_4,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_4(kk,jj)=sum(reshaped_interp_cont_pos_4(kk,:,jj));
    end
end
%%
%+7/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_5(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_5(yy,bb+64.5)=xe_cont_pos_5(bb+64.5)*(3.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_5(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_5(yy,bb+64.5))^2)-((xf_pos_5(yy,bb+64.5))^2)));
        time_diag_cont_pos_5(yy,bb+64.5)=diag_dist_cont_pos_5(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_5(yy,bb+64.5)=time_diag_cont_pos_5(yy,bb+64.5)-time_diag_cont_pos_5(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_5(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_5(dd,:);
end

for hh=1:128
    temp_cont_pos_5=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_5(:,hh),'linear',0);
    reshaped_interp_cont_pos_5(:,hh,:)=reshape(temp_cont_pos_5,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_5(kk,jj)=sum(reshaped_interp_cont_pos_5(kk,:,jj));
    end
end
%%
%+5/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_6(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_6(yy,bb+64.5)=xe_cont_pos_6(bb+64.5)*(2.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_6(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_6(yy,bb+64.5))^2)-((xf_pos_6(yy,bb+64.5))^2)));
        time_diag_cont_pos_6(yy,bb+64.5)=diag_dist_cont_pos_6(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_6(yy,bb+64.5)=time_diag_cont_pos_6(yy,bb+64.5)-time_diag_cont_pos_6(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_6(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_6(dd,:);
end

for hh=1:128
    temp_cont_pos_6=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_6(:,hh),'linear',0);
    reshaped_interp_cont_pos_6(:,hh,:)=reshape(temp_cont_pos_6,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_6(kk,jj)=sum(reshaped_interp_cont_pos_6(kk,:,jj));
    end
end
%%
%+3/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_7(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_7(yy,bb+64.5)=xe_cont_pos_7(bb+64.5)*(1.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_7(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_7(yy,bb+64.5))^2)-((xf_pos_7(yy,bb+64.5))^2)));
        time_diag_cont_pos_7(yy,bb+64.5)=diag_dist_cont_pos_7(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_7(yy,bb+64.5)=time_diag_cont_pos_7(yy,bb+64.5)-time_diag_cont_pos_7(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_7(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_7(dd,:);
end

for hh=1:128
    temp_cont_pos_7=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_7(:,hh),'linear',0);
    reshaped_interp_cont_pos_7(:,hh,:)=reshape(temp_cont_pos_7,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_7(kk,jj)=sum(reshaped_interp_cont_pos_7(kk,:,jj));
    end
end
%%
%+1/2
for yy=1:length(zf)
 for bb=-63.5:1:63.5
        xe_cont_pos_8(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_pos_8(yy,bb+64.5)=xe_cont_pos_8(bb+64.5)*(.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_pos_8(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_pos_8(yy,bb+64.5))^2)-((xf_pos_8(yy,bb+64.5))^2)));
        time_diag_cont_pos_8(yy,bb+64.5)=diag_dist_cont_pos_8(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_pos_8(yy,bb+64.5)=time_diag_cont_pos_8(yy,bb+64.5)-time_diag_cont_pos_8(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_pos_8(dd,:)=timeArray2(dd,:)+time_delay_cont_pos_8(dd,:);
end

for hh=1:128
    temp_cont_pos_8=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_pos_8(:,hh),'linear',0);
    reshaped_interp_cont_pos_8(:,hh,:)=reshape(temp_cont_pos_8,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_pos_8(kk,jj)=sum(reshaped_interp_cont_pos_8(kk,:,jj));
    end
end
%%
%-15/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_1(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_1(yy,bb+64.5)=xe_cont_neg_1(bb+64.5)*(-7.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_1(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_1(yy,bb+64.5))^2)-((xf_neg_1(yy,bb+64.5))^2)));
        time_diag_cont_neg_1(yy,bb+64.5)=diag_dist_cont_neg_1(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_1(yy,bb+64.5)=time_diag_cont_neg_1(yy,bb+64.5)-time_diag_cont_neg_1(yy,65);
    end

end
for dd=1:length(timeArray)
    time_withDelays_cont_neg_1(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_1(dd,:);
end

for hh=1:128
    temp_cont_neg_1=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_1(:,hh),'linear',0);
    reshaped_interp_cont_neg_1(:,hh,:)=reshape(temp_cont_neg_1,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_1(kk,jj)=sum(reshaped_interp_cont_neg_1(kk,:,jj));
    end
end
%%
%-13/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_2(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_2(yy,bb+64.5)=xe_cont_neg_2(bb+64.5)*(-6.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_2(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_2(yy,bb+64.5))^2)-((xf_neg_2(yy,bb+64.5))^2)));
        time_diag_cont_neg_2(yy,bb+64.5)=diag_dist_cont_neg_2(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_2(yy,bb+64.5)=time_diag_cont_neg_2(yy,bb+64.5)-time_diag_cont_neg_2(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_2(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_2(dd,:);
end

for hh=1:128
    temp_cont_neg_2=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_2(:,hh),'linear',0);
    reshaped_interp_cont_neg_2(:,hh,:)=reshape(temp_cont_neg_2,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_2(kk,jj)=sum(reshaped_interp_cont_neg_2(kk,:,jj));
    end
end
%%
%-11/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_3(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_3(yy,bb+64.5)=xe_cont_neg_3(bb+64.5)*(-5.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_3(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_3(yy,bb+64.5))^2)-((xf_neg_3(yy,bb+64.5))^2)));
        time_diag_cont_neg_3(yy,bb+64.5)=diag_dist_cont_neg_3(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_3(yy,bb+64.5)=time_diag_cont_neg_3(yy,bb+64.5)-time_diag_cont_neg_3(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_3(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_3(dd,:);
end

for hh=1:128
    temp_cont_neg_3=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_3(:,hh),'linear',0);
    reshaped_interp_cont_neg_3(:,hh,:)=reshape(temp_cont_neg_3,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_3(kk,jj)=sum(reshaped_interp_cont_neg_3(kk,:,jj));
    end
end
%%
%-9/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_4(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_4(yy,bb+64.5)=xe_cont_neg_4(bb+64.5)*(-4.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_4(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_4(yy,bb+64.5))^2)-((xf_neg_4(yy,bb+64.5))^2)));
        time_diag_cont_neg_4(yy,bb+64.5)=diag_dist_cont_neg_4(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_4(yy,bb+64.5)=time_diag_cont_neg_4(yy,bb+64.5)-time_diag_cont_neg_4(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_4(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_4(dd,:);
end

for hh=1:128
    temp_cont_neg_4=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_4(:,hh),'linear',0);
    reshaped_interp_cont_neg_4(:,hh,:)=reshape(temp_cont_neg_4,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_4(kk,jj)=sum(reshaped_interp_cont_neg_4(kk,:,jj));
    end
end
%%
%-7/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_5(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_5(yy,bb+64.5)=xe_cont_neg_5(bb+64.5)*(-3.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_5(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_5(yy,bb+64.5))^2)-((xf_neg_5(yy,bb+64.5))^2)));
        time_diag_cont_neg_5(yy,bb+64.5)=diag_dist_cont_neg_5(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_5(yy,bb+64.5)=time_diag_cont_neg_5(yy,bb+64.5)-time_diag_cont_neg_5(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_5(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_5(dd,:);
end

for hh=1:128
    temp_cont_neg_5=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_5(:,hh),'linear',0);
    reshaped_interp_cont_neg_5(:,hh,:)=reshape(temp_cont_neg_5,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_5(kk,jj)=sum(reshaped_interp_cont_neg_5(kk,:,jj));
    end
end
%%
%-5/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_6(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_6(yy,bb+64.5)=xe_cont_neg_6(bb+64.5)*(-2.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_6(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_6(yy,bb+64.5))^2)-((xf_neg_6(yy,bb+64.5))^2)));
        time_diag_cont_neg_6(yy,bb+64.5)=diag_dist_cont_neg_6(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_6(yy,bb+64.5)=time_diag_cont_neg_6(yy,bb+64.5)-time_diag_cont_neg_6(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_6(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_6(dd,:);
end

for hh=1:128
    temp_cont_neg_6=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_6(:,hh),'linear',0);
    reshaped_interp_cont_neg_6(:,hh,:)=reshape(temp_cont_neg_6,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_6(kk,jj)=sum(reshaped_interp_cont_neg_6(kk,:,jj));
    end
end
%%
%-3/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_7(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_7(yy,bb+64.5)=xe_cont_neg_7(bb+64.5)*(-1.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_7(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_7(yy,bb+64.5))^2)-((xf_neg_7(yy,bb+64.5))^2)));
        time_diag_cont_neg_7(yy,bb+64.5)=diag_dist_cont_neg_7(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_7(yy,bb+64.5)=time_diag_cont_neg_7(yy,bb+64.5)-time_diag_cont_neg_7(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_7(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_7(dd,:);
end

for hh=1:128
    temp_cont_neg_7=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_7(:,hh),'linear',0);
    reshaped_interp_cont_neg_7(:,hh,:)=reshape(temp_cont_neg_7,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_7(kk,jj)=sum(reshaped_interp_cont_neg_7(kk,:,jj));
    end
end
%%
%-1/2
for yy=1:length(zf)

for bb=-63.5:1:63.5
        xe_cont_neg_8(yy,bb+64.5)=((anecoicCyst.elementSpacingMM)/1000)*bb;
        xf_neg_8(yy,bb+64.5)=xe_cont_neg_8(bb+64.5)*(-.5)*((anecoicCyst.elementSpacingMM)/1000);
        diag_dist_cont_neg_8(yy,bb+64.5)=sqrt(zf(yy)^2 + (((xe_cont_neg_8(yy,bb+64.5))^2)-((xf_neg_8(yy,bb+64.5))^2)));
        time_diag_cont_neg_8(yy,bb+64.5)=diag_dist_cont_neg_8(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont_neg_8(yy,bb+64.5)=time_diag_cont_neg_8(yy,bb+64.5)-time_diag_cont_neg_8(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont_neg_8(dd,:)=timeArray2(dd,:)+time_delay_cont_neg_8(dd,:);
end

for hh=1:128
    temp_cont_neg_8=interp1(timeArray2(:,hh),squeeze(anecoicCyst_data(:,hh,1:16:128)),time_withDelays_cont_neg_8(:,hh),'linear',0);
    reshaped_interp_cont_neg_8(:,hh,:)=reshape(temp_cont_neg_8,[2353,1,8]);
end 

for jj=1:8
    for kk=1:2353
            zone_interp_cont_16_neg_8(kk,jj)=sum(reshaped_interp_cont_neg_8(kk,:,jj));
    end
end
%%
zone_interp_cont_16=zeros(2353,128);

for mm=0:7
    zone_interp_cont_16(:,16*mm+1)=zone_interp_cont_16_neg_1(:,mm+1);
    zone_interp_cont_16(:,16*mm+2)=zone_interp_cont_16_neg_2(:,mm+1);
    zone_interp_cont_16(:,16*mm+3)=zone_interp_cont_16_neg_3(:,mm+1);
    zone_interp_cont_16(:,16*mm+4)=zone_interp_cont_16_neg_4(:,mm+1);
    zone_interp_cont_16(:,16*mm+5)=zone_interp_cont_16_neg_5(:,mm+1);
    zone_interp_cont_16(:,16*mm+6)=zone_interp_cont_16_neg_6(:,mm+1);
    zone_interp_cont_16(:,16*mm+7)=zone_interp_cont_16_neg_7(:,mm+1);
    zone_interp_cont_16(:,16*mm+8)=zone_interp_cont_16_neg_8(:,mm+1);
    zone_interp_cont_16(:,16*mm+9)=zone_interp_cont_16_pos_8(:,mm+1);
    zone_interp_cont_16(:,16*mm+10)=zone_interp_cont_16_pos_7(:,mm+1);
    zone_interp_cont_16(:,16*mm+11)=zone_interp_cont_16_pos_6(:,mm+1);
    zone_interp_cont_16(:,16*mm+12)=zone_interp_cont_16_pos_5(:,mm+1);
    zone_interp_cont_16(:,16*mm+13)=zone_interp_cont_16_pos_4(:,mm+1);
    zone_interp_cont_16(:,16*mm+14)=zone_interp_cont_16_pos_3(:,mm+1);
    zone_interp_cont_16(:,16*mm+15)=zone_interp_cont_16_pos_2(:,mm+1);
    zone_interp_cont_16(:,16*mm+16)=zone_interp_cont_16_pos_1(:,mm+1);
end 

figure;
imagesc(20*log10(abs(hilbert(zone_interp_cont_16))),[30 80])
title('16 parallel receive beams')
colormap('gray')
