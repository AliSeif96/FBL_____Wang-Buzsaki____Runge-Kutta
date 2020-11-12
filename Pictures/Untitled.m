clc
clear
A=dlmread('temp.txt');


x=A(:,1);
y1=A(:,2);
y2=A(:,3);
y3=A(:,4);




plot(x,y1)
hold on
plot(x,y2)
hold on
plot(x,y3)
%axis([30 50 -70 30]);
