%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gerchberg and Saxton Algorithm for retrieveing the phase of the desired
%intensity pattern and projecting onto the spatial light modulator
%author - Hisay Lama
%email ID - hisaylama@gmail.com
%The part of the code, that includes the projection of image onto the
%spatial light modulator is borrowed from the 'ots1.0.1' written by Volpe
%et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,  clc,
close all,
%Input beam 
M = 800; %side length (pixel unit)
N = 600; % side length (pixel unit)
psize = 20e-6; %pxiel size (m)
x=([0.5:1:M-0.5] - M/2)*psize;
y =([0.5:1:N-0.5] - N/2)*psize;
[X,Y]=meshgrid(x,y);
w=250*psize; % Radius of the incident beam (px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input beam projected on diffracted optical element
%u_in = circ((sqrt((X-M/4).^2 + (Y-N/4).^2))./w);
u_in = circ((sqrt(X.^2 + Y.^2))/w); %Ampliture of the input beam
I_in = abs(u_in.^2);
figure(1), %irradiance image
imagesc(x,y,I_in);
xlabel('x (m)'); ylabel('y (m)');
colormap('gray');
%axis square;
axis xy;
%I = double(I_in);
PH = 2.*pi.*(rand([N,M])); %Generation of the random phase
%u_in = u_in.*exp(-1i.*PH); %Resultant input beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the target image with only amplitude i.e. without phase
S_L = input('Enter the number of spots =  ');
P = input('Enter the inter-spot separation in pixel =  ');
w = 1*psize; %radius of the spot (in pixel unit) 
u_target = Multiple_Spot(M,N,S_L,P*psize);
x=([0.5:1:M-0.5] - M/2)*psize;
y =([0.5:1:N-0.5] - N/2)*psize;
[X,Y]=meshgrid(x,y);
I_target = abs(u_target.^2);
figure,
imagesc(x,y,I_target);
%axis off;
axis xy;
colormap('gray'); xlabel('x (px)'); ylabel('y (px)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting of the GS algorithm iteration for geenrating phase mask
f = input('focal length = '); %focal length in pixel units
z = 0; %distance from the focus
lambda = 1064e-9; %input('lambda = '); %light wavelength [m]
k = 2.*pi/lambda; %wave vector
[fx,fy] = meshgrid((1/psize).*[-M/2:1:M/2-1]/M,(1/psize)*[-N/2:1:N/2-1]/N);%frequency size
x1 = lambda*f*fx;
y1 = lambda*f*fy;
Phase_hot = PH;
N_max = input('maximum number of iteration = ');

for n = 1:N_max %number of iteration
    Efocus = exp(1i*k*(2*f+z + x1.^2 + y1.^2))/(1i*lambda*f)...
        .*fftshift(fft2(abs(u_in).*exp(-1j*Phase_hot).*exp((-1i*pi*z/(lambda*f^2))*(X.^2+Y.^2))));
    %Efocus = fftshift(fft2(abs(u_in).*exp(-1j*Phase_hot).*exp((-1i*pi*z/(lambda*f^2))*(X.^2+Y.^2))));
    phase = mod(angle(Efocus), 2*pi);     
    E_in = abs(u_target).*exp(-1j*(phase));
    E_inv = ifftshift(ifft2(ifftshift(E_in)));
    Phase_hot = mod(angle(E_inv), 2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the phase mask
%phase_z = 2*pi*(x.^2 + y.^2).*z/(lambda*f^2);
Phase = mod(Phase_hot, 2*pi);
%Phase = Phase_hot;
figure,
imshow(mat2gray((Phase)+ 1.2*pi))%))%%
%title('Phase Plot')
axis equal tight off
box off
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the resulaant image after introducing the 
R = fftshift(ifft2(fftshift(exp(-1j*(angle(E_inv))))));
figure; imshow(mat2gray(abs(R)));
colorbar()
%imshow(mat2gray(abs(I1)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%Proojecting onto the spatial light modulator. 
%%spatial light modulator SLM and the screen SCRNUM.
%
scrnum = 2; %Screen Number
ScreenPos = get(0,'MonitorPositions');
ScreenWidth = ScreenPos(scrnum,3)-ScreenPos(scrnum,1)+1;
ScreenHeight = ScreenPos(scrnum,4)-ScreenPos(scrnum,2)+1;
%
%The size of the phase mask and of the screen should be the same       
%
if [M N]~=[ScreenHeight ScreenWidth]
   warning('The size of the phase mask and of the screen should be the same.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%slmfig = hologram display figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slmfig = figure( ...
         'Name', 'SLM Screen', ...
         'NumberTitle', 'off', ...
         'MenuBar', 'none', ...
         'ToolBar', 'none', ...
         'Units', 'pixels', ...
         'Resize', 'off', ...
         'Position', [ScreenPos(scrnum,1) ScreenPos(scrnum,2) ScreenWidth ScreenHeight] ...
          );
axes('Units','Normalized','Position',[0,0,1,1],'Visible','off')
figure(slmfig)
cla %clear the axes
image(Phase_hot/(2*pi)*255)
axis equal tight off
colormap(gray(255))
drawnow()
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%closes all figures but keeps hologram.          
fh = findall(0,'type','figure').';
close(fh(fh~=slmfig));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%Uing the Matlab graphics
%MATLAB defines the figure Position property as a vector.
%[left bottom width height]
%MATLAB does not measure the window border when placing the figure; the Position property defines only the internal active area of the figure window.
%The figure's Units property determines the units of the values used to specify the position on the screen. Possible values for the Units property are.
%set(gcf,'Units')
%[ inches | centimeters | normalized | points | {pixels} | 
%characters]
%with pixels being the default.
%
%Determining Screen Size
%get(0,'ScreenSize')
%ans =
%   1 1 1152 900
%In this case, the screen is 1152 by 900 pixels. MATLAB returns the ScreenSize in the units determined by the root Units property. For example,
%set(0,'Units','normalized')
%normalizes the values returned by ScreenSize:
%get(0,'ScreenSize')
%ans =
%    0 0 1 1
%Defining the figure Position in terms of the ScreenSize in normalized units makes the specification independent of variations in screen size. This is useful if you are writing an M-file that is to be used on different computer systems.