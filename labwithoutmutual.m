clc;
clear all;
disp( ' formation of bus admittance matrix ')
nb = input ( 'enter the no of buses' );
nl = input ( 'enter the no of lines' );
for j= 1:nl
sb(j)=input ('enter the starting bus');
eb(j)=input ('enter the ending bus ');
z(j) =input ('enter the input values');
b(j) =input ('enter the values');
end
y=zeros(nl,nb);
for i=1:nl
k=sb(i);
m=eb(i);
y(i,k)=1;
y(i,m)=-1;
end
u=zeros(nl,nl);
for i=1:nl
u(i,i)=1/z(i)+b(i)/2;
end
t=y';
x=t*u*y;
x
y
u
t