%%code for creating the multiple circular spots in the 2D plane in the maltab
 %author - Hisay Lama
 %email - hisaylama@gmail.com
 %spacing betwwen the dots = L/10
 %I = irradiance
 %size of the lattice S_L(i,j) = (4,4)
 %size of the frame = [L/2, L/2 - dx]
function [I1] = Multiple_Spot(M, N, S_L, P)
psize = 20e-6;
w=1*psize; %half width of the spot (px) 
x1=([0.5:1:M-0.5] - M/2)*psize;
y1 =([0.5:1:N-0.5] - N/2)*psize;
[X1,Y1]=meshgrid(x1,y1);
I1 = 0;
for i = 1:1:S_L
    I = 0;
    for j = 1:1:S_L
        %u1=circ((sqrt((X1- M/2 - i.*P).^2 + (Y1 - N/2 - j.*P).^2))./w);
        u1=circ((sqrt((X1 - i.*P).^2 + (Y1 - j.*P).^2))./w);
        I = I + abs(u1.^2); 
    end
    I1= I1 + I;  
end


end
 