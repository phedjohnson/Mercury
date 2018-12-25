%Philip Johnson
%Simple plotter for simulation results.
%As written, it grabs about half of 
%the output files for a given field
%variable from ReorganizeData.cpp

%The ReorganizeData.cpp program must
%be executed before this plotter program
%will work properly

close all
clc
clear all

%Choices for the plotted field:
%Ma = Mach number
%sos = speed of sound
%vx = x velocity component
%vy = y velocity component
%rho = density
%p = pressure

field = 'vx'
fname0 = sprintf('%s_t0.csv',field)
fname2 = sprintf('%s_t2.csv',field)
fname4 = sprintf('%s_t4.csv',field)
fname6 = sprintf('%s_t6.csv',field)
fname8 = sprintf('%s_t8.csv',field)
fname10 = sprintf('%s_t10.csv',field)

tname0 = sprintf('%s (t0)',field)
tname2 = sprintf('%s (t2)',field)
tname4 = sprintf('%s (t4)',field)
tname6 = sprintf('%s (t6)',field)
tname8 = sprintf('%s (t8)',field)
tname10 = sprintf('%s (t10)',field)


size = csvread('size_params.csv');
Mx = size(1,1)
My = size(1,2)

%GEO file contains x,y coordinates
%of cell centroids.
GEO = csvread('XY.csv');
for i = 1:Mx
    for j = 1:My
            x(i,j) = GEO(My*(i-1) + j, 1);
            y(i,j) = GEO(My*(i-1) + j, 2);
    end
end
clear GEO;

DATA = csvread(fname0);
data_min = min(DATA);
data_max = max(DATA);
for i = 1:Mx
    for j = 1:My
        data_0(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;
DATA = csvread(fname2);
data_min = min(data_min,min(DATA));
data_max = max(data_max,max(DATA));
for i = 1:Mx
    for j = 1:My
        data_2(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;
DATA = csvread(fname4);
data_min = min(data_min,min(DATA));
data_max = max(data_max,max(DATA));
for i = 1:Mx
    for j = 1:My
        data_4(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;

DATA = csvread(fname6);
data_min = min(data_min,min(DATA));
data_max = max(data_max,max(DATA));
for i = 1:Mx
    for j = 1:My
        data_6(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;

DATA = csvread(fname8);
data_min = min(data_min,min(DATA));
data_max = max(data_max,max(DATA));
for i = 1:Mx
    for j = 1:My
        data_8(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;

DATA = csvread(fname10);
data_min = min(data_min,min(DATA));
data_max = max(data_max,max(DATA));
for i = 1:Mx
    for j = 1:My
        data_10(i,j) =DATA(My*(i-1) + j, 1);
    end
end
clear DATA;

figure(100)
surface(x,y,data_0)
title(tname0)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(102)
surface(x,y,data_2)
title(tname2)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(104)
surface(x,y,data_4)
title(tname4)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(106)
surface(x,y,data_6)
title(tname6)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(108)
surface(x,y,data_8)
title(tname8)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(110)
surface(x,y,data_10)
title(tname10)
shading flat
caxis([data_min data_max])
colormap jet
colorbar

figure(200)
contourf(x,y,data_0,15)
axis([0 2 0 1])
title(tname0)
colorbar

figure(202)
contourf(x,y,data_2,15)
axis([0 2 0 1])
title(tname2)
colorbar

figure(204)
contourf(x,y,data_4,15)
axis([0 2 0 1])
title(tname4)
colorbar

figure(206)
contourf(x,y,data_6,15)
axis([0 2 0 1])
title(tname6)
colorbar

figure(208)
contourf(x,y,data_8,15)
axis([0 2 0 1])
title(tname8)
colorbar

figure(210)
contourf(x,y,data_10,15)
axis([0 2 0 1])
title(tname10)
colorbar
