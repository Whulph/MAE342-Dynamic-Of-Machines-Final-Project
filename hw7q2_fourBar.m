close all
clc
clear

omega = -1; %rad/s

a = 15;
b = 50;
c = 41.5;
d = 38.7922672707;
h = 3.06;
lineAQ = 1.65;
delta = deg2rad(25);
gamma = deg2rad(31);
x0 = [7*pi/4;3*pi/2;0;0;0;0];

m2 = 1;
m3 = 5;
m4 = 2.33;

IG2G2 = (1/12)*m2*a^2;
IG3G3 = 7;
IG4G4 = (1/12)*m4*c^2;
%Establish time vector
tMax = 2*pi/abs(omega); %max time
nSamples = 200; %number of theta2 angles to input
if omega < 0
    theta2 = linspace(2*pi,0,nSamples+1);
else
    theta2 = linspace(0,2*pi,nSamples+1);
end

%In this case, the initial and final angles are the same (0 and 2pi).
%Remove the last index to ensure smooth motion as it loops.
t = linspace(0,tMax,nSamples+1);

%Remove the last index to ensure smooth motion as it loops.
theta2(end) = [];
t(end) = [];

%Set the tolerance for the magnitude of imaginary number that means the
%solution does not exist.
imagTol = 1e-6;

nSamples = length(t);


omega2 = omega + zeros(1,nSamples);
alpha2 = zeros(1,nSamples);

theta3 = NaN(size(theta2));
omega3 = NaN(size(theta2));
alpha3 = NaN(size(theta2));

theta4 = NaN(size(theta2));
omega4 = NaN(size(theta2));
alpha4 = NaN(size(theta2));

FO2x = NaN(size(theta2));
FO2y = NaN(size(theta2));
FAx = NaN(size(theta2));
FAy = NaN(size(theta2));
FBx = NaN(size(theta2));
FBy = NaN(size(theta2));
FO4x = NaN(size(theta2));
FO4y = NaN(size(theta2));
Tin = NaN(size(theta2));

%Calculate the joint velocities over time. We can do so numerically
fTol = 1e-12;
xTol = 1e-12;
maxIt = 1000;
for i=1:nSamples
    %x(1) = theta3;
    %x(2) = d;
    %x(3) = omega3;
    %x(4) = dd
    %x(5) = alpha3
    %x(6) = ddd
    th2 = theta2(i);
    om2 = omega2(i);
    al2 = alpha2(i);
    f = @(x) [  a*cos(th2) + b*cos(x(1)) - c*cos(x(2)) - d*cos(delta);...
        a*sin(th2) + b*sin(x(1)) - c*sin(x(2)) - d*sin(delta);...
        a*om2*sin(th2) + b*x(3)*sin(x(1)) - c*x(4)*sin(x(2));...
        a*om2*cos(th2) + b*x(3)*cos(x(1)) - c*x(4)*cos(x(2));...
        -a*al2*sin(th2) - a*om2^2*cos(th2) - b*x(5).*sin(x(1)) - b*x(3).^2.*cos(x(1)) + c*x(6).*sin(x(2)) + c*x(4).^2.*cos(x(2));...
        a*al2*cos(th2) - a*om2^2*sin(th2) + b*x(5).*cos(x(1)) - b*x(3).^2.*sin(x(1)) - c*x(6).*cos(x(2)) + c*x(4).^2.*sin(x(2))];
    
    dx = 1e-6;
    dx1 = [dx;0;0;0;0;0];
    dx2 = [0;dx;0;0;0;0];
    dx3 = [0;0;dx;0;0;0];
    dx4 = [0;0;0;dx;0;0];
    dx5 = [0;0;0;0;dx;0];
    dx6 = [0;0;0;0;0;dx];
    
    df = @(x) [ (f(x + dx1) - f(x - dx1))/(2*dx),...
        (f(x + dx2) - f(x - dx2))/(2*dx),...
        (f(x + dx3) - f(x - dx3))/(2*dx),...
        (f(x + dx4) - f(x - dx4))/(2*dx),...
        (f(x + dx5) - f(x - dx5))/(2*dx),...
        (f(x + dx6) - f(x - dx6))/(2*dx)];
    
    if i == 1
        %x0 = [0;pi/2;0;omega2(i);0;0];
    else
        x0 = [theta3(i-1);theta4(i-1);omega3(i-1);omega4(i-1);alpha3(i-1);alpha4(i-1)];
    end
    
    xf = newtonRaphsonND( f,df,x0,fTol,xTol,maxIt,[]);
    theta3(i) = xf(1);
    try theta4(i) = xf(2);
    catch
        error('newtonRaphsonND returned NaN in loop %i',i)
    end
    omega3(i) = xf(3);
    omega4(i) = xf(4);
    alpha3(i) = xf(5);
    alpha4(i) = xf(6);
    
    
end

%Use results to animate crank-slider motion.
%Calculate the location of each joint over time.
%Within a for-loop, we will update the XData and YData fields of each plot
%object, then save that frame to a video file, effectively turning the
%figure into a video.

RO4complex = d*exp(1i*delta);
RO4 = [real(RO4complex);imag(RO4complex)];

RAcomplex = a*exp(1i*theta2);
RA = [real(RAcomplex);imag(RAcomplex)];

RBcomplex = RAcomplex + b*exp(1i*theta3);
RB = [real(RBcomplex);imag(RBcomplex)];

RO2 = [0;0];

RPcomplex = RAcomplex + h*exp(1i*(theta3 + gamma));
RP = [real(RPcomplex);imag(RPcomplex)];

VAcomplex = 1i*omega2.*a.*exp(1i*theta2);
AAcomplex = (1i*omega2).^2.*a.*exp(1i*theta2) + 1i*alpha2.*a.*exp(1i*theta2);

VPcomplex = VAcomplex + 1i*omega3.*h.*exp(1i*(theta3 + gamma));
VP = [real(VPcomplex);imag(VPcomplex)];

APcomplex = AAcomplex + (1i*omega3).^2.*h.*exp(1i*(theta3 + gamma)) + 1i*alpha3.*h.*exp(1i*(theta3 + gamma));
AP = [real(APcomplex);imag(APcomplex)];


RQcomplex = RAcomplex + lineAQ*exp(1j*(theta3+gamma));
RQ = [real(RQcomplex);imag(RQcomplex)];

VQcomplex = VAcomplex + 1j*omega3.*lineAQ.*exp(1j*(theta3+gamma));
VQ = [real(VQcomplex);imag(VQcomplex)];

AQcomplex = AAcomplex - lineAQ*omega3.^2.*exp(1j*(theta3+gamma));
AQ = [real(AQcomplex);imag(AQcomplex)];

RG2complex = (a/2)*exp(1j*theta2);
RG2 = [real(RG2complex);imag(RG2complex)];

VG2complex = 1i*omega2.*(a/2).*exp(1i*theta2);
AG2complex = (1i*omega2).^2.*(a/2).*exp(1i*theta2) + 1i*alpha2.*(a/2).*exp(1i*theta2);
AG2 = [real(AG2complex);imag(AG2complex)];

RG4complex = RO4complex + (c/2)*exp(1j*theta4);
VG4complex = 1j*omega4.*(c/2).*exp(1j*theta4);
AG4complex = (1j*omega4).^2.*(c/2).*exp(1j*theta4) + 1j*alpha4.*(c/2).*exp(1j*theta4);
AG4 = [real(AG4complex);imag(AG4complex)];

VPMax = max(sqrt(sum(VP.^2)));

APMax = max(sqrt(sum(AP.^2)));

h3 = figure;
subplot(2,2,1);
plot(theta2,real(VPcomplex),'linewidth',2)
hold on
grid on
xlim([0,2*pi])
xticks(0:pi/2:2*pi)
xlabel('\theta_2 (rad)')
ylim([-VPMax,VPMax])
ylabel('V_{P_x}, length/sec^2')
title('x component of V_P')
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

subplot(2,2,2)
plot(theta2,imag(VPcomplex),'linewidth',2)
hold on
grid on
xlim([0,2*pi])
xticks(0:pi/2:2*pi)
xlabel('\theta_2 (rad)')
ylim([-VPMax,VPMax])
ylabel('V_{P_y}, length/sec^2')
title('y component of V_P')
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

subplot(2,2,3);
plot(theta2,real(APcomplex),'linewidth',2)
hold on
grid on
xlim([0,2*pi])
xticks(0:pi/2:2*pi)
xlabel('\theta_2 (rad)')
ylim([-APMax,APMax])
ylabel('A_{P_x}, length/sec^2')
title('x component of A_P')
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

subplot(2,2,4)
plot(theta2,imag(APcomplex),'linewidth',2)
hold on
grid on
xlim([0,2*pi])
xticks(0:pi/2:2*pi)
xlabel('\theta_2 (rad)')
ylim([-APMax,APMax])
ylabel('A_{P_y}, length/sec^2')
title('y component of A_P')
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})


h3.Position(2) = 50;
h3.Position(3) = 700;
h3.Position(4) = 700;


for i=1:nSamples
    
    k31 = (a/2)*sin(theta2(i));
    k32 = (a/2)*cos(theta2(i));
    k33 = (a/2)*sin(theta2(i));
    k34 = (a/2)*cos(theta2(i));
    k43 = lineAQ*sin(theta3(i)+gamma);
    k44 = lineAQ*cos(theta3(i)+gamma);
    k45 = b*sin(theta3(i)) - lineAQ*sin(theta3(i)+gamma);
    k46 = c*cos(theta3(i)) - lineAQ*cos(theta3(i)+gamma);
    k85 = -(c/2)*sin(theta4(i));
    k86 = -(c/2)*cos(theta4(i));
    k87 = -(c/2)*sin(theta4(i));
    k88 = -(c/2)*cos(theta4(i));
    
    A = [1 0 -1 0 0 0 0 0 0;
        0 -1 0 1 0 0 0 0 0;
        k31 k32 k33 k34 0 0 0 0 -1;
        0 0 1 0 -1 0 0 0 0;
        0 0 0 -1 0 1 0 0 0;
        0 0 k43 k44 k45 k46 0 0 0;
        0 0 0 0 1 0 -1 0 0;
        0 0 0 0 0 -1 0 1 0;
        0 0 0 0 k85 k86 k87 k88 0;
        ];
    
    B = [
        m2*AG2(1,i);
        m2*AG2(2,i);
        IG2G2*alpha2(i);
        m3*AQ(1,i);
        m3*AQ(2,i);
        IG3G3*alpha3(i);
        m4*AG4(1,i);
        m4*AG4(2,i);
        IG4G4*alpha4(i);
        ];
    
    %C = [FO2x;FO2y;FAx;FAy;FBx;FBy;FO4x;FO4y;Tin];
    
    C = linsolve(A,B);
    
    FO2x(i) = C(1);
    FO2y(i) = C(2);
    FAx(i) = C(3);
    FAy(i) = C(4);
    FBx(i) = C(5);
    FBy(i) = C(6);
    FO4x(i) = C(7);
    FO4y(i) = C(8);
    Tin(i) = C(9);
    
end

figure
hold on 
subplot(1,3,1)
plot(t,Tin)
grid on
xlabel('Time')
ylabel('Input Torque (Tin)')

subplot(1,3,2)
plot(theta2,Tin)
grid on
xlabel('\theta_2')
ylabel('Input Torque (Tin)')

subplot(1,3,3)
plot(omega2,Tin)
grid on
xlabel('\omega_2')
ylabel('Input Torque (Tin)')
