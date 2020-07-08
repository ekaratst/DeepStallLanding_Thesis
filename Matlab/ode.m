
%-------------------------------------AOA----------------------------------------------------------------------------------
 %-------------------------------------ODE----------------------------------------------------------------------------------
% x_degrees = linspace(-10,180);
% index_figure = 1;
% for i = 1:1
%     x_degrees = [0, 5, 10, 15];
%     x = x_degrees(i)*pi/180; %aoa
%     t_interval = [0,10]; 
%     aircraft_velocity = 8;
%     pitch_rate = 0;
%     height = 0;
%     horizontal_distance = 0;
%     for j = 1:1
%         angle_of_descent = x_degrees(j)*pi/180;
%         pitch_angle = x + angle_of_descent;
%         sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));
%         CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
%         CLreg = CL0 + CLalpha .* x;
%         CLdsl = 2*sign(x).*(sin(x).^2) .* cos(x);
%         CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
%         CDreg = CD0 + K .* (CL0 + CLalpha .* x) .^2;
%         CDdsl = 2*sign(x).*(sin(x).^3);
%         CM = CM0 + CMalpha * x;
%         %time interval
%         t_interval = [0,10];
%         %solution
%         [t,y] = ode45(@(t,Y) odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x) , t_interval , [aircraft_velocity, angle_of_descent, pitch_angle, pitch_rate, height, horizontal_distance]);
%         
%         figure(index_figure)
%         subplot(3,2,1);
%         plot(y(:,6),y(:,5)); %vetical distance(y-axis)-horizontal distance(x-axis)
%         xlabel('Horizontal Distance [m]') 
%         ylabel('Vetical Distance [m]') 
%         title(sprintf('Simulated DSL Trajectory(AOA=%d, Angle of descent=%d, Pitch angle=%d)',x_degrees(i), x_degrees(j), x_degrees(i)+x_degrees(j) ))
%         legend("8 m/s")
%         index_figure = index_figure + 1;
%     end
%     
% end

%-------------------------------------AOA-----------------------------------------------------------------------------------
% plot(t,y(:,1),'b',t,y(:,2),'r');

% figure(1)
% subplot(2,1,1);
% plot(x_degrees_aoa, CL_aoa);
% % plot(x_degrees_aoa, sin_test);
% hold on
% plot(x_degrees_aoa, CLreg_aoa);
% plot(x_degrees_aoa, CLdsl_aoa);
% %xlim([-10 90])
% ylim([-2 2])
% title('C_L - AOA')
% xlabel('Angel of Attack [deg]') 
% ylabel('Lift Coefficient') 
% legend({'CL','CL_r_e_g','CL_D_S_L'},'Location','northeast')
% hold off
% 
% subplot(2,1,2);
% plot(x_degrees_aoa, CD_aoa);
% hold on
% plot(x_degrees_aoa, CDreg_aoa);
% plot(x_degrees_aoa, CDdsl_aoa)
% xlim([-10 90])
% %ylim([-5 10])
% title('C_D - AOA')
% xlabel('Angel of Attack [deg]') 
% ylabel('Drag Coefficient') 
% legend({'CD','CD_r_e_g','CD_D_S_L'},'Location','northeast')
% hold off

%------------------------------------ODE------------------------------------------------------------------------------------
figure_number = 1;
aoa = 0;
gramma = 0;
delta = [0,-10,-20,-30, -50, -70];
for j = delta
    figure(figure_number)
    subplot(3,2,1);
    hold on
    for i = [8,12,16,20]
        [t,y] = solve_ode(i,j,aoa,gramma);
        plot(y(:,6),y(:,5)); %vetical distance(y-axis)-horizontal distance(x-axis)
    end
    xlabel('Horizontal Distance [m]') 
    ylabel('Vetical Distance [m]') 
    title('Simulated DSL Trajectory')
    legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
    hold off

    subplot(3,2,2);
    hold on
    for i = [8,12,16,20]
        [t,y] = solve_ode(i,j,aoa,gramma);
        plot(t,y(:,1))
    end
    xlabel('Time') 
    ylabel('V [m/s]') 
    title('DSL Velocity')
    legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
    hold off

    subplot(3,2,3);
    hold on
    for i = [8,12,16,20]
        [t,y] = solve_ode(i,j,aoa,gramma);
        plot(t,y(:,2))
    end
    xlabel('Time') 
    ylabel('gramma') 
    title('DSL Angle of Descent')
    legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
    hold off

    subplot(3,2,4);
    hold on
    for i = [8,12,16,20]
        [t,y] = solve_ode(i,j,aoa,gramma);
        plot(t,y(:,3))
    end
    xlabel('Time') 
    ylabel('theta') 
    title('DSL Pitch Angle')
    legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
    hold off

    subplot(3,2,5);
    hold on
    for i = [8,12,16,20]
        [t,y] = solve_ode(i,j,aoa,gramma);
        plot(t,y(:,4))
    end
    xlabel('Time') 
    ylabel('q') 
    title('DSL Pitch rate')
    legend('8 m/s', '12 m/s', '16 m/s', '20 m/s')
    hold off

%     sgtitle(sprintf('AOA=%d, Angle of descent=%d, Pitch angle=%d',0,0,0))
    sgtitle(sprintf('Elevator Angle=%d \n AOA=%d, Angle of descent=%d, Pitch angle=%d',j, aoa, gramma, aoa+gramma))
    figure_number = figure_number + 1;
end

function dYdt = odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x)
dYdt = [    (T*cos(x)) - (1/2 * density_air * S .* Y(1).^2 .* (CD + ((CDq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CDdelta.*delta)) - (m*g*sin(Y(2))); 
            (T*sin(x)) + (1/2 * density_air * S .* Y(1).^2 .* (CL + ((CLq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CLdelta.*delta)) - (m*g*cos(Y(2))); 
            (1/2 * density_air * S .* Y(1).^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CMdelta.*delta)) ./ Iyy; 
            Y(3); 
            Y(1)*sin(Y(2)); 
            Y(1)*cos(Y(2))];
end

function [t,y] = solve_ode(acv,elevator_angle,aoa,gramma)
    %-----------------------parameters-------------------------------------------------------------------
    Iyy = 0.388;
    m = 1.5;
    g = 9.81;
    density_air = 1.225;
    S = 0.28;
    mean_chord = 0.31;
    CL0 = 0.062;
    CLalpha = 6.098;
    CD0 = 0.098;
    K = 0.012;
    CM0 = 0.028;
    CMalpha = -0.031;
    CLq = 0;
    CLdelta = -1.72;
    CDq = 0;
    CDdelta = -0.814;
    CMq = -13.1;
    CMdelta = -0.325;
    alpha0 = 20*pi/180;
    M = 50;
    delta = elevator_angle*pi/180;
    T = 0;
    x_degrees = aoa; %<<<-----------------------
    x = x_degrees*pi/180; %aoa
    aircraft_velocity = acv;
    pitch_rate = 0;
    height = 0;
    horizontal_distance = 0;
    angle_of_descent = gramma*pi/180; %<<<-----------------------
    pitch_angle = x + angle_of_descent;
    %-----------------------equations-------------------------------------------------------------------
    sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));
    CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
    CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
    CM = CM0 + CMalpha * x;
    %time interval
    t_interval = [0,10];
    %solution
    [t,y] = ode45(@(t,Y) odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x) , t_interval , [aircraft_velocity, angle_of_descent, pitch_angle, pitch_rate, height, horizontal_distance]);
end     


