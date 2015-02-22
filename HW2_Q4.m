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
    time_withDelays_cont(dd,:)=timeArray(dd,:)+time_delay_cont(dd,:);
end


%for hh=1:128
    delay_interp(:,1,1)=interp1(time_withDelays(:,1),pointTarget.data(:,1,1),timeArray,'linear');
    figure;
    imagesc(delay_interp(:,1,1))
    colormap('gray')
%end 

