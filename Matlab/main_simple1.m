t_span = [0 12];      %Time period
Y_init = [  0 
            1 
            1];     %Initial condition
[T1,Y1] = ode45(@simple1,t_span,Y_init);
plot(T1,Y1(:,1),'-', T1,Y1(:,2),'-.', T1,Y1(:,3),':');