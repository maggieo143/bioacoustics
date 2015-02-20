zone_size=round(max(size(pointTarget.data))/5);

zone1=pointTarget.data(1:zone_size,:,:);
zone2=pointTarget.data(zone_size:2*zone_size,:,:);
zone3=pointTarget.data(2*zone_size:3*zone_size,:,:);
zone4=pointTarget.data(3*zone_size:4*zone_size,:,:);
zone5=pointTarget.data(4*zone_size:end,:,:);

center=(zone_size/2)+1;

center_zone1=center;
center_zone2=center+zone_size;
center_zone3=center+2*zone_size;
center_zone4=center+3*zone_size;
center_zone5=center+4*zone_size;

0:1/(pointTarget.samplingRateMHz*(10^6)):(size(pointTarget.data,1)-1)/(pointTarget.samplingRateMHz*(10^6));

for zf=[(time(center_zone1)*1540)/2 (time(center_zone2)*1540)/2 (time(center_zone3)*1540)/2 (time(center_zone4)*1540)/2 (time(center_zone5)*1540)/2]
    for bb=-63.5:1:63.5
        xe(bb+64.5)=((pointTarget.elementSpacingMM)/1000)*bb;
        diag_dist(bb+64.5)=sqrt(zf^2 + (xe(bb+64.5))^2);
        time_diag(bb+64.5)=diag_dist(bb+64.5)/1540; 
        time_delay(bb+64.5)=time_diag(bb+64.5)-time_diag(65)
    end
end 

