%%gauss seidel main program
clc;
clear all;
close all;
%         |  From |  To   |   R   |   X   |   B/2 |
%         |  Bus  | Bus   |       |       |       |
%linedata = [ 1      2       0.08    0.24     0;
            % 1         3       0    0.06   0;
            % 2      3       0.06    0.018    0 ];
         
%ybus=ybusgs(linedata) %create function Ybus with the ybus code
disp( ' formation of bus admittance matrix ') 
nb = input ( 'enter the no of buses' ); 
nl = input ( 'enter the no of lines' ); 
for j= 1:nl 
    sb(j)=input ('enter the starting bus'); 
    eb(j)=input ('enter the ending bus  '); 
    z(j) =input ('enter the input values'); 
    b(j) =input ('enter the  values'); 
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
ybus=x;
 
% %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
busdata = [ 1     2    1.05     0     0.0    0    0      0      0      0;
            2     3    0        0     400    250    0.5      0.25      0   0;
            3     3    0        0     0.0    0    0    0.3   0   0 ];
bus = busdata(:,1);         % Bus number.
type = busdata(:,2);        % Type of Bus 1-Slack, 2-PV, 3-PQ.
V = busdata(:,3);           % Initial Bus Voltages.
th = busdata(:,4);          % Initial Bus Voltage Angles.
GenMW = busdata(:,5);       % PGi, Real Power injected into the buses.
GenMVAR = busdata(:,6);     % QGi, Reactive Power injected into the buses.
LoadMW = busdata(:,7);      % PLi, Real Power Drawn from the buses.
LoadMVAR = busdata(:,8);    % QLi, Reactive Power Drawn from the buses.
Qmin = busdata(:,9);        % Minimum Reactive Power Limit
Qmax = busdata(:,10);       % Maximum Reactive Power Limit
nbus = max(bus);            % To get no. of buses
P = GenMW - LoadMW;         % Pi = PGi - PLi, Real Power at the buses.
Q = GenMVAR - LoadMVAR;     % Qi = QGi - QLi, Reactive Power at the buses.
Vprev = V;
toler = 1;                  % Tolerence.
iteration = 1;              % iteration starting
while (toler > 0.00001)     % Start of while loop
    for i = 2:nbus
        sumyv = 0;
        for k = 1:nbus
            if i ~= k
                sumyv = sumyv + ybus(i,k)* V(k);  % Vk * Yik
            end
        end
        if type(i) == 2     % Computing Qi for PV bus
            Q(i) = -imag(conj(V(i))*(sumyv + ybus(i,i)*V(i)));
            if (Q(i) > Qmax(i)) || (Q(i) < Qmin(i))  % Checking for Qi Violation.
                if Q(i) < Qmin(i)   % Whether violated the lower limit.
                    Q(i) = Qmin(i);
                else    % No, violated the upper limit.
                    Q(i) = Qmax(i);
                end
                type(i) = 3;  % If Violated, change PV bus to PQ bus.
            end
        end
        V(i) = (1/ybus(i,i))*((P(i)-j*Q(i))/conj(V(i)) - sumyv); % Compute Bus Voltages.
        if type(i) == 2 % For PV Buses, Voltage Magnitude remains same, but Angle changes.           
            V(i) =abs(Vprev(i))*cos(angle(V(i)))+j*abs(Vprev(i))*sin(angle(V(i)));
        end
    end
    iteration = iteration + 1;      % Increment iteration count.
    toler = max(abs(abs(V) - abs(Vprev)));     % Calculate tolerance.
    Vprev = V; % Vprev is required for next iteration,  V(i) = pol2rect(abs(Vprev(i)), angle(V(i)));
end     % End of while loop / Iteration
 
iteration;       % Total iterations.
V               % Bus Voltages in Complex form.
Vmag = abs(V)   % Final Bus Voltages.
Ang = 180/pi*angle(V)    % Final Bus Voltage Angles in Degree.

