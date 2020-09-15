function dy = simple1(t,y)
dy = zeros(3,1);
y1 = y(1);
y2 = y(2);
y3 = y(3);
dy(1) = y2*y3;
dy(2) = -y1*y3;
dy(3) = -0.51*y1*y2;
end
