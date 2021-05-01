clc, close all

%link lengths
a = 38;
b = 41.5;
c = 39.3;
d = 40.1;
e = 55.8;
f = 39.4;
g = 36.7;
h = 65.7;
p = 49;
o = 50;
k = 61.9;
l = 7.8;
m = 15;
n = sqrt(a^2+l^2);

triL1 = 49.53;
triL2 = 18.795;

delta1 = deg2rad(86.268);
delta2 = deg2rad(99.096);

phi1 = deg2rad(56.73);
phi2 = deg2rad(39.99);




%masses FROM BALSA WOOD
mm = 1.2;
mo = 3.63;
mk = 4.46;
mc = 2.89;
mf = 2.9;
mb = 14.45;
mg = 13.56;

%define FN
FN = 1*(mm + mo + mk + mc + mf + mb + mg);

%calculate moments of inertia
IGmGm = (1/12)*mm*m^2;
IGoGo = (1/12)*mo*o^2;
IGkGk = (1/12)*mk*k^2;
IGcGc = (1/12)*mc*c^2;
IGfGf = (1/12)*mf*f^2;
IGbGb = 9637.85;
IGgGg = 8087.17;


%Point A
RAcomplex = m*exp(1j*theta22);
RA = [real(RAcomplex);imag(RAcomplex)];

VAcomplex = 1j*omega22.*m.*exp(1j*theta22);
VA = [real(VAcomplex);imag(VAcomplex)];

AAcomplex =  (-omega22.^2 + 1j*alpha22).*m.*exp(1j*theta22);
AA = [real(AAcomplex);imag(AAcomplex)];

%COM of link m
Rcommcomplex =  m/2*exp(1j*theta22);
Rcomm = [real(Rcommcomplex);imag(Rcommcomplex)];

Vcommcomplex =  1j*omega22.*(m/2).*exp(1j*theta22);
Vcomm = [real(Vcommcomplex);imag(Vcommcomplex)];

Acommcomplex =  (-omega22.^2 + 1j*alpha22).*(m/2).*exp(1j*theta22);
Acomm = [real(Acommcomplex);imag(Acommcomplex)];

%Point B
RBcomplex = RA + o*exp(1j*theta33);
RB = [real(RBcomplex);imag(RBcomplex)];

VBcomplex = VA + 1j*omega33.*o.*exp(1j*theta33);
VB = [real(VBcomplex);imag(VBcomplex)];

ABcomplex = AA + (-omega33.^2 + 1j*alpha33).*o.*exp(1j*theta33);
AB = [real(ABcomplex);imag(ABcomplex)];

%COM of link 0
Rcom0complex = RA + o/2*exp(1j*theta33);
Rcom0 = [real(Rcom0complex);imag(Rcom0complex)];

Vcom0complex = VA + 1j*omega33.*(o/2).*exp(1j*theta33);
Vcom0 = [real(Vcom0complex);imag(Vcom0complex)];

Acom0complex = AA + (-omega33.^2 + 1j*alpha33).*(o/2).*exp(1j*theta33);
Acom0 = [real(Acom0complex);imag(Acom0complex)];

%Point 04
R04complex = -n*exp(1j*theta1);
R04 = [real(R04complex);imag(R04complex)];

%COM of link b
Rcombcomplex = R04 + b/2*exp(1j*theta44);
Rcomb = [real(Rcombcomplex);imag(Rcombcomplex)];

Vcombcomplex = 1j*omega44.*(b/2).*exp(1j*theta44);
Vcomb = [real(Vcombcomplex);imag(Vcombcomplex)];

Acombcomplex = (-omega44.^3 + 1j*alpha44).*(b/2).*exp(1j*theta44);
Acomb = [real(Acombcomplex);imag(Acombcomplex)];

%Point C
RCcomplex = RA + k*exp(1j*theta55);
RC = [real(RCcomplex);imag(RCcomplex)];

VCcomplex = VA + 1j*omega55.*k.*exp(1j*theta55);
VC = [real(VCcomplex);imag(VCcomplex)];

ACcomplex = AA + (-omega55.^2 + 1j*alpha55).*k.*exp(1j*theta55);
AC = [real(ACcomplex);imag(ACcomplex)];

%COM of link k
Rcomkcomplex =  RA + k/2*exp(1j*theta55);
Rcomk = [real(Rcomkcomplex);imag(Rcomkcomplex)];

Vcomkcomplex =  VA + 1j*omega55.*(k/2).*exp(1j*theta55);
Vcomk = [real(Vcomkcomplex);imag(Vcomkcomplex)];

Acomkcomplex =  AA + (-omega55.^2 + 1j*alpha55).*(k/2).*exp(1j*theta55);
Acomk = [real(Acomkcomplex);imag(Acomkcomplex)];

%COM of link c
Rcomccomplex = R04 + c*exp(1j*theta66);
Rcomc = [real(Rcomccomplex);imag(Rcomccomplex)];

Vcomccomplex = 1j*omega66.*c.*exp(1j*theta66);
Vcomc = [real(Vcomccomplex);imag(Vcomccomplex)];

Acomccomplex = (-omega66.^2 + 1j*alpha66).*c.*exp(1j*theta66);
Acomc = [real(Acomccomplex);imag(Acomccomplex)];

%Point E
REcomplex = RC + g*exp(1j*theta88);
RE = [real(REcomplex);imag(REcomplex)];

VEcomplex = VC + 1j*omega88.*g.*exp(1j*theta88);
VE = [real(VEcomplex);imag(VEcomplex)];

AEcomplex = AC + (-omega88.^2 + 1j*alpha88).*g.*exp(1j*theta88);
AE = [real(AEcomplex);imag(AEcomplex)];

%COM of link g
Rcomgcomplex = RC + g/2*exp(1j*theta88);
Rcomg = [real(Rcomgcomplex);imag(Rcomgcomplex)];

Vcomgcomplex = VC + 1j*omega88.*(g/2).*exp(1j*theta88);
Vcomg = [real(Vcomgcomplex);imag(Vcomgcomplex)];

Acomgcomplex = AC + (-omega88.^2 + 1j*alpha88).*(g/2).*exp(1j*theta88);
Acomg = [real(Acomgcomplex);imag(Acomgcomplex)];


%Point D
RDcomplex = R04 + d*exp(1j*(theta44-delta1));
RD = [real(RDcomplex);imag(RDcomplex)];

VDcomplex = 1j*omega44.*d.*exp(1j*(theta44-delta1));
VD = [real(VDcomplex);imag(VDcomplex)];

ADcomplex = (-omega44.^2 + 1j*alpha44).*d.*exp(1j*(theta44-delta1));
AD = [real(ADcomplex);imag(ADcomplex)];

%COM of link d
Rcomdcomplex = R04 + d/2*exp(1j*(theta44-delta1));
Rcomd = [real(Rcomdcomplex);imag(Rcomdcomplex)];

Vcomdcomplex = 1j*omega44.*(d/2).*exp(1j*(theta44-delta1));
Vcomd = [real(Vcomdcomplex);imag(Vcomdcomplex)];

Acomdcomplex = (-omega44.^2 + 1j*alpha44).*(d/2).*exp(1j*(theta44-delta1));
Acomd = [real(Acomdcomplex);imag(Acomdcomplex)];

%COM of link f
Rcomfcomplex = RD + f/2*exp(1j*theta77);
Rcomf = [real(Rcomfcomplex);imag(Rcomfcomplex)];

Vcomfcomplex = VD + 1j*omega77.*(f/2).*exp(1j*theta77);
Vcomf = [real(Vcomfcomplex);imag(Vcomfcomplex)];

Acomfcomplex = AD + (-omega77.^2 + 1j*alpha77).*(f/2).*exp(1j*theta77);
Acomf = [real(Acomfcomplex);imag(Acomfcomplex)];

%Point P
RPcomplex = RC + p*exp(1j*(2*pi-(delta2-theta88)));
RP = [real(RPcomplex);imag(RPcomplex)];

VPcomplex = VC + 1j*(-omega88).*p.*exp(1j*(2*pi-(delta2-theta88)));
VP = [real(VPcomplex);imag(VPcomplex)];

APcomplex = AC + (omega88.^2 + 1j*alpha88).*p.*exp(1j*(2*pi-(delta2-theta88)));
AP = [real(APcomplex);imag(APcomplex)];

subplot(3,1,1)
plot(RP(1,:),RP(2,:))
xlabel('X position (length)')
ylabel('Y position (length)')
title('Foot Position Path')
grid on

subplot(3,1,3)
hold on
plot(t,VP(1,:))
plot(t,VP(2,:))
xlabel('time (sec)')
% xticks(0:pi/4:2*pi)
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('velocity (unit/s)')
legend('VPx','VPy')
title('Foot velocities over time')
grid on

subplot(3,1,2)
hold on
plot(t,RP(1,:))
plot(t,RP(2,:))
xlabel('time (sec)')
% xticks(0:pi/4:2*pi)
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('Position (unit)')
legend('RPx','RPy')
title('Foot position over time')
grid on




%{
h = figure;

hold on

linkm = plot([R02(1,1),RA(1,1)],[R02(2,1),RA(2,1)],'linewidth',2);
linko = plot([RA(1,1),RB(1,1)],[RA(2,1),RB(2,1)],'linewidth',2);
linkb = plot([R04(1,1),RB(1,1)],[R04(2,1),RB(2,1)],'linewidth',2);
linkk = plot([RA(1,1),RC(1,1)],[RA(2,1),RC(2,1)],'linewidth',2);
linkc = plot([R04(1,1),RC(1,1)],[R04(2,1),RC(2,1)],'linewidth',2);
linkg = plot([RC(1,1),RE(1,1)],[RC(2,1),RE(2,1)],'linewidth',2);
linkf = plot([RD(1,1),RE(1,1)],[RD(2,1),RE(2,1)],'linewidth',2);
linkd = plot([R04(1,1),RD(1,1)],[R04(2,1),RD(2,1)],'linewidth',2);


O2 = plot(R02(1,1),R02(2,1),'ko','linewidth',2);
O4 = plot(R04(1,1),R04(2,1),'ko','linewidth',2);

linkn = plot([R02(1,1),R04(1,1)],[R02(2,1),R04(2,1)],'k--');

tt = title(sprintf('t = %3.3f',t(1)));

%Resize the window to contain every possible configuration
xlim([R02(1,1)-2*max(a,c),R04(1,1)+2*max(a,c)])
ylim([R02(2,1)-2*max(a,c),R02(2,1)+2*max(a,c)])
axis equal
ax = gca;
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';

drawnow

%Initialize the videowriter object. This is a structure that will store our
%video frames and then export them as a .avi video.
v = VideoWriter('4barEx1');
v.FrameRate = nSamples/tMax;

open(v);

frame = getframe(h);
writeVideo(v,frame);

%For each time step, update the x and y position of each joint, then save
%the figure to the video structure.
for i=2:nSamples
    
    linkm.XData = [R02(1,i),RA(1,i)];
    linkm.YData = [R02(2,i),RA(2,i)];
    
    linko.XData = [RA(1,i),RB(1,i)];
    linko.YData = [RA(2,i),RB(2,i)];
    
%     linkb.XData = [R04(1,i),RB(1,i)];
%     linkb.YData = [R04(2,i),RB(2,i)];

    linkk.XData = [RA(1,i),RC(1,i)];
    linkk.YData = [RA(2,i),RC(2,i)];
    
%     linkc.XData = [R04(1,i),RC(1,i)];
%     linkc.YData = [R04(2,i),RC(2,i)];
    
    linkg.XData = [RC(1,i),RE(1,i)];
    linkg.YData = [RC(2,i),RE(2,i)];
    
    linkf.XData = [RD(1,i),RE(1,i)];
    linkf.YData = [RD(2,i),RE(2,i)];
    
%     linkd.XData = [R04(1,i),RD(1,i)];
%     linkd.YData = [R04(2,i),RD(2,i)];
    
    tt = title(sprintf('t = %3.3f',t(i)));

    
    drawnow
    
    frame = getframe(h);
    writeVideo(v,frame);
end

close(v);
%}

for i=1:nSamples
    
    k31 = (m/2)*sin(theta22(i));
    k32 = -(m/2)*cos(theta22(i));
    k33 = -(m/2)*sin(theta22(i));
    k34 = (m/2)*cos(theta22(i));
    k65 = (o/2)*cos(theta33(i));
    k66 = -(o/2)*sin(theta33(i));
    k67 = (o/2)*sin(theta33(i));
    k68 = -(o/2)*cos(theta33(i));
    k99 = (k/2)*sin(theta55(i));
    k910 = -(k/2)*cos(theta55(i));
    k911 = -(k/2)*sin(theta55(i));
    k912 = (k/2)*cos(theta55(i));
    k1213 = -(c/2)*sin(theta66(i));
    k1214 = (c/2)*cos(theta66(i));
    k1215 = (c/2)*sin(theta66(i));
    k1216 = -(c/2)*cos(theta66(i));
    k1517 = -(f/2)*sin(theta77(i));
    k1518 = (f/2)*cos(theta77(i));
    k1519 = -(f/2)*sin(theta77(i));
    k1520 = (f/2)*cos(theta77(i));
    k185 = -triL1*sin(theta44(i)-delta1+phi1) + b*sin(theta44(i));
    k186 = triL1*cos(theta44(i)-delta1+phi1) - b*cos(theta44(i));
    k1819 = triL1*sin(theta44(i)-delta1+phi1) - d*sin(theta44(i)-delta1);
    k1820 = -triL1*cos(theta44(i)-delta1+phi1) + d*cos(theta44(i)-delta1);
    k1821 = triL1*sin(theta44(i)-delta1+phi1);
    k1822 = -triL1*cos(theta44(i)-delta1+phi1);
    k2117 = -triL2* + g*sin(theta88(i));
    k2118 = triL2*cos(theta88(i)-delta2+phi2) - g*cos(theta88(i));
    k2123 = triL2*sin(theta88(i)-delta2+phi2);
    k2124 = -triL2*cos(theta88(i)-delta2+phi2);
    kfn = -triL2*cos(theta88(i)-delta2+phi2) + p*cos(theta88(i)-delta2);
    
    A = [1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        k31 k32 k33 k34 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 k65 k66 k67 k68 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 k99 k910 k911 k912 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 k1213 k1214 k1215 k1216 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k1517 k1518 k1519 k1520 0 0 0 0 0 0 0;
        0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0;
        0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0;
        0 0 0 0 k185 k186 0 0 0 0 0 0 0 0 0 0 0 0 k1819 k1820 k1821 k1822 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k2117 k2118 0 0 0 0 k2123 k2124 0 0 0;
        0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0;
        ];
    
    B = [
        mm*Acomm(1,i);
        mm*Acomm(2,i);
        IGmGm*alpha22(i);
        mo*Acom0(1,i);
        mo*Acom0(2,i);
        IGoGo*alpha33(i);
        mk*Acomk(1,i);
        mk*Acomk(2,i);
        IGkGk*alpha55(i);
        mc*Acomc(1,i);
        mc*Acomc(2,i);
        IGcGc*alpha66(i);
        mf*Acomf(1,i);
        mf*Acomf(2,i);
        IGfGf*alpha77(i);
        mb*Acomb(1,i);
        mb*Acomb(2,i);
        IGbGb*alpha44(i);
        mg*Acomg(1,i);
        mg*Acomg(2,i)-FN;
        IGgGg*alpha88(i)-(kfn*FN);
        0;
        0;
        0;
        0;
        0;
        0;
        ];
    %{
    C = [1-FO2x;
         2-FO2y;
         3-FAxm;
         4-FAym;
         5-FBx;
         6-FBy;
         7-FAxo;
         8-FAyo;
         9-FAxk;
         10-FAyk;
         11-FCxk;
         12-FCyk;
         13-FCxc;
         14-FCyc
         15-FO4xc;
         16-FO4yc;
         17-FEx;
         18-FEy;
         19-FDx;
         20-FDy;
         21-FO4xb;
         22-FO4yb;
         23-FCxg;
         24-FCyg;
         25-F04xn;
         26-F04yn;
         27-Tin];
        %}
        C = linsolve(A,B);
        
        FO2x(i) = C(1);
        FO2y(i) = C(2);
        FAxm(i) = C(3);
        FAym(i) = C(4);
        FBx(i) = C(5);
        FBy(i) = C(6);
        FAxo(i) = C(7);
        FAyo(i) = C(8);
        FAxk(i) = C(9);
        FAyk(i) = C(10);
        FCxk(i) = C(11);
        FCyk(i) = C(12);
        FCxc(i) = C(13);
        FCyc(i) = C(14);
        FO4xc(i) = C(15);
        FO4yc(i) = C(16);
        FEx(i) = C(17);
        FEy(i) = C(18);
        FDx(i) = C(19);
        FDy(i) = C(20);
        FO4xb(i) = C(21);
        FO4yb(i) = C(22);
        FCxg(i) = C(23);
        FCyg(i) = C(24);
        FO4xn(i) = C(25);
        FO4yn(i) = C(26);
        Tin(i) = C(27);
        
end

figure
hold on
plot(t,FO2x/12)
plot(t,FAxm/12)
plot(t,FBx/12)
plot(t,FAxo/12)
plot(t,FAxk/12)
plot(t,FCxk/12)
plot(t,FCxc/12)
plot(t,FO4xc/12)
plot(t,FEx/12)
plot(t,FDx/12)
plot(t,FO4xb/12)
plot(t,FCxg/12)
plot(t,FO4xn/12)
xlabel('time (sec)')
ylabel('Force (lbf)')
legend('FO2x','FAxm','FBx','FAxo','FAxk','FCxk','FCxc','FO4xc','FEx','FDx','F04xb','FCxg', 'F04xn')
title('Solution to x direction forces over time')
grid on

figure
hold on
plot(t,FO2y/12)
plot(t,FAym/12)
plot(t,FBy/12)
plot(t,FAyo/12)
plot(t,FAyk/12)
plot(t,FCyk/12)
plot(t,FCyc/12)
plot(t,FO4yc/12)
plot(t,FEy/12)
plot(t,FDy/12)
plot(t,FO4yb/12)
plot(t,FCyg/12)
plot(t,FO4yn/12)
xlabel('time (sec)')
ylabel('Force (lbf)')
legend('FO2y','FAym','FBy','FAyo','FAyk','FCyk','FCyc','FO4yc','FEy','FDy','F04yb','FCyg', 'F04yn')
title('Solution to y direction forces over time')
grid on

figure
subplot(3,1,1)
hold on
plot(t,Tin/24)
xlabel('time (sec)')
ylabel('Torque (lbf/ft)')
title('Torque Over Time')
grid on

subplot(3,1,2)
hold on
plot(theta22,Tin/24)
xlabel('Angular Crank Displacement (rad)')
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('Torque (lbf/ft)')
title('Input Torque Versus Angular Crank Displacement')
grid on

subplot(3,1,3)
hold on
plot(omega22,Tin/24)
xlabel('Angular Crank Velocity (rad/sec)')
ylabel('Torque (lbf/ft)')
title('Input Torque Versus Angular Crank Velocity')
grid on




