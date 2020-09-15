convert_to_degree = 180/pi;
t_span = [0 10];      %Time period
Y_init = [   0       % Velocity(V)
             0       % Angle of descent(gramma) -30*pi/180
             0       % Pitch rate(q)
             0       % Pitch angle(theta) 20*pi/180
             0       % height(h)
             0];     % horizontal distance(r)

figure(1)
for i = 1:4
    vel = [8,12,16,20];
    Y_init(1) = vel(i);
    [T1,Y1] = ode45(@longitudinal_equation,t_span,Y_init);  
    plot(Y1(:,6),Y1(:,5)); 
    hold on 
end
xlabel('Horizontal Distance [m]') 
ylabel('Vetical Distance [m]') 
title("Simulated DSL Trajectory)")
legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
hold off

figure(2)
for i = 1:4
    vel = [8,12,16,20];
    Y_init(1) = vel(i);
    [T1,Y1] = ode45(@longitudinal_equation,t_span,Y_init);  
    plot(T1,Y1(:,1)); 
    hold on 
end
xlabel('Time[sec]') 
ylabel('V[m/s]') 
title("Simulated DSL Velocity")
legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
hold off

figure(3)
for i = 1:4
    vel = [8,12,16,20];
    Y_init(1) = vel(i);
    [T1,Y1] = ode45(@longitudinal_equation,t_span,Y_init);  
    plot(T1,Y1(:,2)*convert_to_degree); 
    hold on 
end
xlabel('Time[sec]') 
ylabel('\gamma') 
title("Simulated DSL Angle of descent(\gamma)")
legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
hold off


