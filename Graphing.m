clc
close all

theta66 

h1 = figure;
subplot(3,1,1)
hold on
plot(theta22,theta33)
plot(theta22,theta33)
plot(theta22,theta55)
plot(theta22,theta66)
plot(theta22,theta77)
plot(theta22,theta88)
xlabel('crank angle (radians)')
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('absolute angle (rad)')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8')
title('Newton-Raphson solution to angles over time')
grid on


subplot(3,1,2)
hold on
plot(theta22,omega33)
plot(theta22,omega33)
plot(theta22,omega55)
plot(theta22,omega66)
plot(theta22,omega77)
plot(theta22,omega88)
xlabel('crank angle (rad)')
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('angular velocity (rad/s)')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8')
title('Newton-Raphson solution to velocities over time')
grid on


subplot(3,1,3)
hold on
plot(theta22,alpha33)
plot(theta22,alpha33)
plot(theta22,alpha55)
plot(theta22,alpha66)
plot(theta22,alpha77)
plot(theta22,alpha88)
xlabel('crank angle (rad)')
xticks(0:pi/4:2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
ylabel('angular accelration (rad/s^2)')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8')
title('Newton-Raphson solution to accelerations over time')
grid on

h1.Position(3) = 450;
drawnow