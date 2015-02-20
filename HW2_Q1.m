field_init(-1)

for zz=[1 10 20 30 40 50 60 70]
%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a phased array B-mode system scans an image
%
%  This script assumes that the field_init procedure has been called
%
%  Example by Joergen Arendt Jensen, Nov. 28, 1995.
%  Generate the transducer apertures for send and receive

f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength
element_height=5/1000;   %  Height of element [m]
kerf=0.1/1000;           %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]

%  Generate aperture for emission
emit_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);

%  Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

%  Generate aperture for reception
receive_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);

%  Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

% Do phased array imaging
%point_position=[0 0 70]/1000;   %  Position of the point to be imaged 
no_lines=50;                    %  Number of A-lines in image
sector=20 * pi/180;             %  Size of image sector
d_theta=sector/no_lines;        %  Increment in angle for 90 deg. image

%  Pre-allocate some storage
image_data=zeros(800,no_lines);

theta= -sector/2; 
    
for i=1:no_lines
    
   %   Set the focus for this direction
  xdc_focus (emit_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
  xdc_focus (receive_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
  
  %   Calculate the received response

  point_position=[0 0 zz]/1000;   %  Position of the point to be imaged
  [v, t1]=calc_scat(emit_aperture, receive_aperture, point_position, 1); 
  %end
  
  %  Store the result
  image_data(1:max(size(v)),i)=v;
  times(i) = t1;
  
  %  Steer in another angle
  theta = theta + d_theta;
end

%create samples vector
for ii=1:length(times)
samples(ii)=times(ii)*fs;
end 

%preallocated zeros matrix for point spread function envelope
point_spread_envelope=zeros(max(samples)+max(size(image_data)),min(size(image_data)));

%point spread function envelope
for kk=1:min(size(image_data))
    point_spread_envelope(samples(kk):samples(kk)+max(size(image_data))-1,kk)=image_data(1:max(size(image_data)),kk);
end

cLow=min(min(point_spread_envelope))
cHigh=max(max(point_spread_envelope))

%image of point spread function envelope
figure;
imagesc(20*log10(abs(hilbert(point_spread_envelope))),[cLow,cHigh])
colormap('gray')

figure;
imagesc(point_spread_envelope,[cLow,cHigh])
colormap('gray')
end 


% As the the depth of the scatter changes, and becomes further away from
% the focal point, the 