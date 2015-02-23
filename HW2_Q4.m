for zz=1:length(pointTarget.data)
    zf(zz,1)=(time(zz)*1540)/2;
end 

for yy=1:length(zf)
    for bb=-63.5:1:63.5
        xe_cont(yy,bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        diag_dist_cont(yy,bb+64.5)=sqrt(zf(yy)^2 + (xe_cont(yy,bb+64.5))^2);
        time_diag_cont(yy,bb+64.5)=diag_dist_cont(yy,bb+64.5)/1540;  
    end
    
    for bb=-63.5:1:63.5
    time_delay_cont(yy,bb+64.5)=time_diag_cont(yy,bb+64.5)-time_diag_cont(yy,65);
    end
end

for dd=1:length(timeArray)
    time_withDelays_cont(dd,:)=timeArray2(dd,:)+time_delay_cont(dd,:);
end


for hh=1:128
    delay_interp_cont(:,hh,:)=interp1(time_withDelays_cont(:,hh),squeeze(pointTarget.data(:,hh,:)),timeArray2(:,hh),'linear');
end 

cLow3=min(min(min(delay_interp_cont)));
cHigh3=max(max(max(delay_interp_cont)));
figure;
imagesc(20*log10(abs(hilbert(delay_interp_cont(:,:)))),[cLow3,cHigh3])
colormap('gray')
