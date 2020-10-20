convert_to_degree = 180/pi;
elevator_angle = -45;

t_span = [0 15];      %Time period
Y_init = [   0       % Velocity(V)
             0       % Angle of descent(gramma) -30*pi/180
             0       % Pitch rate(q)
             0       % Pitch angle(theta) 20*pi/180
             15       % height(h)
             0];     % horizontal distance(r)

figure(1)
for i = 1:4
    vel = [8,12,16,20];
    Y_init(1) = vel(i);
    [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle),t_span,Y_init);  
    plot(Y1(:,6),Y1(:,5)); 
    hold on 
end
ylim([0 22])
xlabel('Horizontal Distance [m]') 
ylabel('Vetical Distance [m]') 
title("Simulated DSL Trajectory [\delta = -45]")
legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
hold off         
         
figure(2)
for i = 1:4
    vel = [8,12,16,20];
    Y_init(1) = vel(i);
    [T1,Y1] = ode45(@(t,y) adjust_longitudinal_equation(t,y,elevator_angle),t_span,Y_init);  
    plot(Y1(:,6),Y1(:,5)); 
    hold on 
    
end
ylim([0 22])
xlabel('Horizontal Distance [m]') 
ylabel('Vetical Distance [m]') 
title("Simulated DSL Trajectory [\delta = -45]")
legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
hold off

% figure(2)
% for i = 1:4
%     vel = [8,12,16,20];
%     Y_init(1) = vel(i);
%     [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle),t_span,Y_init);  
%     plot(T1,Y1(:,1)); 
%     hold on 
% end
% xlabel('Time[sec]') 
% ylabel('V[m/s]') 
% title("Simulated DSL Velocity [\delta = -45]")
% legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
% hold off
% 
% figure(3)
% for i = 1:4
%     vel = [8,12,16,20];
%     Y_init(1) = vel(i);
%     [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle),t_span,Y_init);  
%     plot(T1,Y1(:,2)*convert_to_degree); 
%     hold on 
% end
% xlabel('Time[sec]') 
% ylabel('\gamma') 
% title("Simulated DSL Angle of descent(\gamma) [\delta = -45]")
% legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
% hold off
% 
% figure(4)
% for i = 1:7
%     elevator_angle = [-15, -25, -35, -45, -55, -65, -75];
%     Y_init(1) = 15;
%     [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle(i)),t_span,Y_init);  
%     plot(Y1(:,6),Y1(:,5)); 
%     hold on 
% end
% ylim([0 20])
% xlabel('Horizontal Distance [m]') 
% ylabel('Vetical Distance [m]') 
% title("Simulated DSL Trajectory [Each Elevator Angle]")
% legend('-15', '-25', '-35', '-45', '-55', '-65', '-75')
% hold off




% figure(5)
% Y_init(1) = 15;
% if T1 >=0
%     elevator_angle = -45;
%     [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle(i)),t_span,Y_init);  
% else
%     elevator_angle = -65;
%     [T1,Y1] = ode45(@(t,y) longitudinal_equation(t,y,elevator_angle(i)),t_span,Y_init);  
% end
% plot(Y1(:,6),Y1(:,5)); 
% 
% elevator_angle = [0, 1, 2.5, 5, 5.2, 10];
% p7vals = [pi/9, 0.3, exp(-1), sqrt(2), -3];  %one fewer
% num_interval = length(p7vals);
% t = cell(num_interval,1);
% x = cell(num_interval, 1);
% for idx = 1 : num_interval
%     p7 = p7vals(idx);
%     [t{idx},x{idx}] = ode45(@(t,x) try_eqns(t, x, p7), elevator_angle(idx:idx+1), x0);
%     x0 = x{idx}(end,:);
% end
% 



