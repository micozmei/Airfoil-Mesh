%NACA 4 digit airfoil CFD
clear all; close all; clc;
axis equal

tic
Airfoil = '5412';
chordlength = 1;
gridpoints = 300;
radial_density = 50;
AngleofAttack = 5;  %degrees

%% Airfoil Drawing
disp('Drawing Airfoil')

m  = str2double(Airfoil(1));
p  = str2double(Airfoil(2));
xx = str2double(Airfoil(3:4));

%define parameters
M = m/100;
P = p/10;
T = xx/100;
c = chordlength;
step = c/gridpoints;

%Mean Camber Line:
x = zeros(gridpoints,1);
for i = 1:gridpoints
    x(i) = 0.5*cos(2*pi*i/gridpoints)+0.5;
end
x = sort(x);

%x = [0:step:c].';
y_c = [0];
dy_c = [0];
for j = 1:gridpoints
   if x(j) >=0 && x(j) < P.*c;
        y_c(j,1)  = M/(P.^2).*x(j).*(2.*P - x(j)/c);
        dy_c(j,1) = 2.*M/(P.^2).*(P - x(j)/c);  
   else x(j) >= P.*c && x(j) <= c;
        y_c(j,1)  = M.*((c-x(j))/((1-P).^2)).*( 1+x(j)/c - 2.*P );
        dy_c(j,1) = 2.*M/((1 - P).^2).*(x(j)/c - P);    
   end
end

%Angle of Attack 
alpha = -AngleofAttack*pi/180;
Rotation = ([cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[x,y_c].').';
x = Rotation(:,1);
y_c = Rotation(:,2);
hold on

%Thinkness Distribution
a_0 = 0.2969;
a_1 = -0.126;
a_2 = -0.3516;
a_3 = 0.2843;
a_4 = -0.1015;   %trailing edge pointy
%a_4 = -0.1036;  %trailing edge closed

y_t = T/0.2.*(  a_0.*x.^(0.5) + a_1.*x + a_2.*x.^2 + a_3.*x.^3 + a_4.*x.^4);

theta = atan(dy_c);
x_u = x - y_t.*sin(theta)-0.1;
y_u = y_c + y_t.*cos(theta);

x_L = x + y_t.*sin(theta)-0.1;
y_L = y_c - y_t.*cos(theta);

%% Mesh Point Generation 
disp('Generating Mesh')

% ID matrix for airfoil
% Boundary = zeros(2*(gridpoints),2);
% for i = 1:2*gridpoints 
%     Boundary(i,1) = i;
%     Boundary(i,2) = i+1;
% end
% Boundary(2*(gridpoints),1) = 2*(gridpoints);
% Boundary(2*(gridpoints),2) = 1;

%Start the point matrix with the airfoil points
Points(1:gridpoints,1) = x_u;
Points(1:gridpoints,2) = y_u;
Points( gridpoints+1:2*(gridpoints) , 1 ) = flipud(x_L);
Points( gridpoints+1:2*(gridpoints) , 2 ) = flipud(y_L);
enis = length(Points);

for i = 1:enis
        for j = 1:radial_density
            Points(2*(gridpoints) + radial_density*i-(radial_density-j),:) = (exp(2.8*j/radial_density)).*Points(i,:); 
        end
end

%Throw away unnecessary points
sqweenis = length(Points);
for i = 1:sqweenis
    if abs(Points(sqweenis-i+1,1)-0.5)> 1.05
        Points(sqweenis-i+1,:) = [];
    elseif abs( Points(sqweenis-i+1,2) ) > 0.55
        Points(sqweenis-i+1,:) = [];
    end
end

DT = delaunayTriangulation(Points);
triplot(DT,'k')
plot(x-0.1,y_c,'b-')
plot(x_u,y_u,'r-')
plot(x_L,y_L,'r-')
title(Airfoil)
axis([-0.5,1.5,-0.5,0.5])

%% Finite Element Analysis
disp('Running Finite Volume Analysis')
%compare to "eppler code, mark drela code(XFOIL)"
toc