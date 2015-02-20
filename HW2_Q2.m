%% Intro
time=[0:1:size(pointTarget.data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));

rate_upsample=((pointTarget.samplingRateMHz)*(10^6))*4;

time_upsample=[0:1/4:size(pointTarget.data,1)-1]*(1/((pointTarget.samplingRateMHz)*(10^6)));

figure;
for aa=1:128
    interpolation(:,:,aa)=interp1(time,pointTarget.data(:,:,aa),time_upsample,'linear');
    imagesc(interpolation(:,:,aa))
    colormap('gray')
end

zf=0.04;   %[m]

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

for ee=1:length(time_delay)
    samples(ee)=round(time_delay(ee)*rate_upsample);
end

%% Part A.
interpolation_delay=zeros(max(size(time_upsample)),128,128);

for jj=1:128
    for ff=1:128
        shifted_interp=max(size(interpolation))-samples(ff);
        interpolation_delay(1:shifted_interp,ff,jj)=interpolation(samples(ff)+1:end,ff,jj);
    end
end

figure;
[cLow]=min(min(min(interpolation_delay)));
[cHigh]=max(max(max(interpolation_delay)));
imagesc(interpolation_delay(:,:),[cLow, cHigh])
colormap('gray')
 

 %% Part B.
% figure;
% imagesc(interpolation_delay(:,65),[cLow, cHigh])
% colormap('gray')
% 
% %  POINT TARGET???


%% Part C. 
    for jj=1:128
    for kk=1:9725
            sum_interp(kk,jj)=sum(interpolation_delay(kk,:,jj));
    end
    end
   
%      for ll=1:9725
%         total_interp(ll,1)=sum(sum_interp(ll,:));
%     end
  
cLow=min(min(min(total_interp)));
cHigh=max(max(max(total_interp)));

figure;
imagesc(20*log10(abs(hilbert( sum_interp))),[cLow,cHigh])
colormap('gray')


%%

%channel_matrices(:,jj,:)=squeeze(veraStrct.data(:,jj,:));
%interpolation(:,jj,:)=interp1(time,channel_matrices(:,jj,:),time_upsample,'linear');

