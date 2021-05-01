close all
clc

omega = .08*pi; %rad/s

%Input the lengths of the links. d is grounded. a is driven by
%the actuator.

a = 15;
b = 50;
c = 41.5;
d = 38.7922672707;
x0 = [deg2rad(30);deg2rad(67.11);0;0;0;0];

height = 7.8;
width = 38;
theta1 = pi - atan(height/width);

%First, we must check for Grashof condition. If the linkage is Grashof,
%then the input angle may cycle around the entire unit circle. If not, then
%we must calculate the toggle points and restrict motion to be between
%them.
%
%Put all link lengths in a vector. The shortest one is S, the longest one
%is L.
Lvec = [a,b,c,d];
[minVal,minInd] = min(Lvec);
[maxVal,maxInd] = max(Lvec);
S = minVal;
L = maxVal;
%Now, to identify P and Q, remove the indices associated with S and L. The
%lengths in these indices are P and Q, in either order.
PQinds = 1:4;
PQinds([minInd,maxInd]) = [];
P = Lvec(PQinds(1));
Q = Lvec(PQinds(2));

%Check the Grashof condition. 
if S+L < P+Q
    %Grashof!
    %Set isGrashof to true, and set the minimum and maximum angles of
    %theta2 to be 0 and 2pi, respectively.
    isGrashof = true;
    isSpecial = false;
else
    %Not Grashof!
    %Set isGrashof to false, and then solve for the toggle points of the
    %mechanism. These toggle points will limit the motion to be animated.
    isGrashof = false;
    isSpecial = false;
    error('Mechanism must be Grashof for this analysis to work')
end
    
%Solve for the toggle points of the
%mechanism. These toggle points will limit the motion to be animated.

%From equation 4.37, we can calculate the toggle points of the
%mechanism. There are two solutions; only one of which will be real.
%Thus, we must test both to find which is real.
minCandidates = NaN(1,2);
maxCandidates = NaN(1,2);
%The first solution is the case in which the b-c segments are
%collapsed. This is a candidate for minimum angle. But if the
%proportions are such that (a + d) < (b + c), then this will also be
%the maximum angle, once the crank first extends and then collapses the
%linkage.
arg1 = (a^2 + d^2 - b^2 - c^2)/(2*a*d) + (b*c)/(a*d);
if abs(arg1) <= 1
    minCandidates(1) = acos(arg1);
    maxCandidates(1) = 2*pi-acos(arg1);
end
%The second solution is the case in which the b-c segments are
%extended. This is a candidate for maximum angle. But if the
%proportions are such that (d-a) < (b-c), then this will also be the
%minimum angle, once the crank first collabes and then re-extends the
%linkage.
arg2 = (a^2 + d^2 - (b + c)^2)/(2*a*d);
if abs(arg2) <= 1
    maxCandidates(2) = acos(arg2);
end

%Check our limits and make sure we have the proper number.
if all(~isnan(minCandidates))
    theta2min = max(minCandidates);
elseif any(~isnan(minCandidates))
    theta2min = minCandidates(~isnan(minCandidates));
else
    theta2min = 0;
end

if all(~isnan(maxCandidates))
    theta2max = min(maxCandidates);
elseif any(~isnan(maxCandidates))
    theta2max = maxCandidates(~isnan(maxCandidates));
else
    theta2max = 2*pi;
end

%Establish time vector
tMax = 2*pi/omega; %max time
nSamples = 200; %number of theta2 angles to input
if (theta2min == 0) && (theta2max == 2*pi)
    theta2 = linspace(theta2min,theta2max,nSamples+1);
    %In this case, the initial and final angles are the same (0 and 2pi).
    %Remove the last index to ensure smooth motion as it loops.
    t = linspace(0,tMax,nSamples+1);
else
    theta2 = [linspace(theta2min,theta2max,nSamples/2),linspace(theta2max,theta2min,nSamples/2)];
    t = linspace(0,tMax,nSamples);
end
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

theta1 = zeros(size(theta2)) + theta1;

%Calculate the joint velocities over time. We can do so numerically
fTol = 1e-12;
xTol = 1e-12;
maxIt = 1000;
for i=1:nSamples
    %x(1) = theta3;
    %x(2) = theta4;
    %x(3) = omega3;
    %x(4) = omega4
    %x(5) = alpha3
    %x(6) = alpha4
    th2 = theta2(i);
    om2 = omega2(i);
    al2 = alpha2(i);
    f = @(x) [  a*cos(th2) + b*cos(x(1)) - c*cos(x(2)) + d*cos(theta1(i));...
                a*sin(th2) + b*sin(x(1)) - c*sin(x(2)) + d*sin(theta1(i));...
                -a*om2*sin(th2) - b*x(3)*sin(x(1)) + c*x(4)*sin(x(2));...
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

theta22 = theta2;
omega22 = omega2;
alpha22 = alpha2;
theta33 = theta3;
theta44 = theta4;
omega33 = omega3;
omega44 = omega4;
alpha33 = alpha3;
alpha44 = alpha4;


%Plot a preliminary figure of each angle over time. These are basic data
%that may help us diagnose problems with our system.
h1 = figure;
subplot(3,1,1)
plot(t,theta2)
hold on
plot(t,theta3)
plot(t,theta4)
xlabel('time (s)')
ylabel('absolute angle (rad)')
legend('\theta_2','\theta_3','\theta_4')
title('Newton-Raphson solution to angles over time')


subplot(3,1,2)
plot(t,omega2)
hold on
plot(t,omega3)
plot(t,omega4)
xlabel('time (s)')
ylabel('angular velocity (rad/s)')
legend('\omega_2','\omega_3','\omega_4')
title('Newton-Raphson solution to angular velocities over time')

subplot(3,1,3)
plot(t,alpha2)
hold on
plot(t,alpha3)
plot(t,alpha4)
xlabel('time (s)')
ylabel('angular acceleration (rad/s^2)')
legend('\alpha_2','\alpha_3','\alpha_4')
title('Newton-Raphson solution to angular accelerations over time')

h1.Position(3) = 450;
drawnow

%Use results to animate 4-bar motion.
%Calculate the location of each joint over time.
%Within a for-loop, we will update the XData and YData fields of each plot
%object, then save that frame to a video file, effectively turning the
%figure into a video.

%link "a"
RAcomplex = a*exp(1j*theta2);
RA = [real(RAcomplex);imag(RAcomplex)];

%link "c"
RBcomplex = RA + b*exp(1i*theta3);
RB = [real(RBcomplex);imag(RBcomplex)];

%anchor points
R02 = zeros(2,length(theta2));

R04complex = -d*exp(1i*theta1);
R04 = [real(R04complex);imag(R04complex)];

%Initialize a figure and plot the initial configuration of the system.
h = figure;
link2 = plot([R02(1,1),RA(1,1)],[R02(2,1),RA(2,1)],'linewidth',2);
hold on
link3 = plot([RA(1,1),RB(1,1)],[RA(2,1),RB(2,1)],'linewidth',2);
link4 = plot([RB(1,1),R04(1,1)],[RB(2,1),R04(2,1)],'linewidth',2);
% AAvec = plot([RcomA(1,1),RcomA(1,1) + AcomA(1,1)],[RcomA(2,1), RcomA(2,1) + AcomA(2,1)],'w');
% ABvec = plot([RcomB(1,1),RcomB(1,1) + AcomB(1,1)],[RcomB(2,1), RcomB(2,1) + AcomB(2,1)],'k');
% ACvec = plot([RcomC(1,1),RcomC(1,1) + AcomC(1,1)],[RcomC(2,1), RcomC(2,1) + AcomC(2,1)],'k');
O2 = plot(R02(1,1),R02(2,1),'ko','linewidth',2);
O4 = plot(R04(1,1),R04(2,1),'ko','linewidth',2);
link1 = plot([R02(1,1),R04(1,1)],[R02(2,1),R04(2,1)],'k--');
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
    
    link2.XData = [R02(1,i),RA(1,i)];
    link2.YData = [R02(2,i),RA(2,i)];
    
    link3.XData = [RA(1,i),RB(1,i)];
    link3.YData = [RA(2,i),RB(2,i)];
    
    link4.XData = [RB(1,i),R04(1,i)];
    link4.YData = [RB(2,i),R04(2,i)];
    
    %Update the time in the title.
    tt = title(sprintf('t = %3.3f',t(i)));

    
    drawnow
    
    frame = getframe(h);
    writeVideo(v,frame);
end

close(v);

mv = omega4./omega2;

h2 = figure;
subplot(2,2,1)
plot(t,mv,'linewidth',2)
xlabel('time (s)')
title('Angular velocity ratio m_v')
grid on
ylim([-10,10])

subplot(2,2,2)
plot(theta2,mv,'linewidth',2)
xlabel('\theta_2 (rad)')
title('Angular velocity ratio m_v')
xlim([theta2min,theta2max])
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
grid on
ylim([-10,10])

ma = 1./mv;

subplot(2,2,3)
plot(t,ma,'linewidth',2)
xlabel('time (s)')
title('Mechanical advantage m_a')
grid on
ylim([-10,10])

subplot(2,2,4)
plot(theta2,ma,'linewidth',2)
xlabel('\theta_2 (rad)')
title('Mechanical advantage ratio m_a')
xlim([theta2min,theta2max])
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
grid on

h2.Position = [319.4000 277 874.4000 420];
ylim([-10,10])
