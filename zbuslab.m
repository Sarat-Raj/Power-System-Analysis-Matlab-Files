clear all
clc
disp('Zbus Building Algorithm')
data=dlmread('data.m');
num_bus=max(max(data(:,1),data(:,2)));
sz=size(data);
size=sz(1,1);
buses_added=1;
bus_status=zeros(1,num_bus+1);
for n=1:size
    if(data(n,1)==0)
        temp=data(n,1);
        data(n,1)=data(n,2);
        data(n,2)=temp;
    end
    if(data(n,1)>data(n,2) && data(n,2)~=0)
        temp=data(n,1);
        data(n,1)=data(n,2);
        data(n,2)=temp;
    end
end
for n=1:size
    if(data(n,1)==1 && data(n,2)==0)
        temp1=data(1,:);
        data(1,:)=data(n,:);
        data(n,:)=temp1;
    end
end
for n=1:size
    for m=1:size
        if(data(m,1)==n)
            bus1=data(m,1);
            bus2=data(m,2);
            p_bus1=0;
            p_bus2=0;
            for k=1:num_bus
                if(bus_status(1,k)== bus1 && bus1~=0)
                    p_bus1=1;
                end
                if(bus_status(1,k)== bus2 && bus2~=0)
                    p_bus2=1;
                end
            end
            if(bus_status(1,buses_added)==0 && bus2==0 && p_bus1==0)
                disp('Adding Z=')
                disp(data(m,3)+ i*data(m,4))
                disp('between buses:')
                disp(bus1)
                disp(bus2)
                disp('This impedance is added between a new bus and reference')
                buses_added=buses_added+1;
                bus_status(1,buses_added-1)=bus1;
                if(bus1==1)
                    Zbus(bus1,bus1)=data(m,3)+i*data(m,4)
                else
                    ssz=length(Zbus);
                    Zbus(ssz+1,ssz+1)=data(m,3)+i*data(m,4)
                end
                disp(' ')
                disp(' ')
            elseif(p_bus1==1 && p_bus2==0 && bus2~=0)
                disp('Adding Z=')
                disp(data(m,3)+ i*data(m,4))
                disp('between buses:')
                disp(bus1)
                disp(bus2)
                disp('This impedance is added between a new bus and an existing bus')
                buses_added=buses_added+1;
                bus_status(1,buses_added-1)=bus2;
                size_zbus=length(Zbus);
                for var=1:size_zbus
                    Zbus(size_zbus+1,var)=Zbus(bus1,var);
                    Zbus(var,size_zbus+1)=Zbus(var,bus1);
                end
                Zbus(size_zbus+1,size_zbus+1)=Zbus(bus1,bus1)+ data(m,3)+ i*data(m,4);
                Zbus
                disp(' ')
                disp(' ')
            elseif(p_bus1==1 && p_bus2==0 && bus2==0)
                disp('Adding Z=')
                disp(data(m,3)+ i*data(m,4))
                disp('between buses:')
                disp(bus1)
                disp(bus2)
                disp('This impedance is added between an existing bus and reference')
                size_zbus=length(Zbus);
                Zbus1=Zbus;
                for var=1:size_zbus
                    Zbus1(size_zbus+1,var)=Zbus1(bus1,var);
                    Zbus1(var,size_zbus+1)=Zbus1(var,bus1);
                end
                Zbus1(size_zbus+1,size_zbus+1)=Zbus1(bus1,bus1)+ data(m,3)+ i*data(m,4);
                for var1=1:size_zbus
                    for var2=1:size_zbus
                        Zbus(var1,var2)=Zbus1(var1,var2)- Zbus1(var1,size_zbus+1)*Zbus1(size_zbus+1,var2)/Zbus1(size_zbus+1,size_zbus+1);
                    end
                end
                Zbus
                disp(' ')
                disp(' ')
            elseif(p_bus1==1 && p_bus2==1 && bus2~=0)
                disp('Adding Z=')
                disp(data(m,3)+ i*data(m,4))
                disp('between buses:')
                disp(bus1)
                disp(bus2)
                disp('This impedance is added between two existing buses')
                size_zbus=length(Zbus);
                Zbus1=Zbus;
                for var=1:size_zbus
                    Zbus1(size_zbus+1,var)=Zbus1(bus1,var)-Zbus1(bus2,var);
                    Zbus1(var,size_zbus+1)=Zbus1(var,bus1)-Zbus1(var,bus2);
                end
                Zbus1(size_zbus+1,size_zbus+1)=Zbus1(bus1,bus1)+Zbus1(bus2,bus2)-2*Zbus1(bus1,bus2) + data(m,3)+ i*data(m,4);
                for var1=1:size_zbus
                    for var2=1:size_zbus
                        Zbus(var1,var2)=Zbus1(var1,var2)- Zbus1(var1,size_zbus+1)*Zbus1(size_zbus+1,var2)/Zbus1(size_zbus+1,size_zbus+1);
                    end
                end
                Zbus
                disp(' ')
                disp(' ')
            end
        end
    end
end
bus_status
for var1=1:num_bus
    for var2=1:num_bus
        if(bus_status(1,var2)==var1)
            zvar1=Zbus(var2,:);
            Zbus(var2,:)=Zbus(var1,:);
            Zbus(var1,:)=zvar1;
            zvar2=Zbus(:,var2);
            Zbus(:,var2)=Zbus(:,var1);
            Zbus(:,var1)=zvar2;
            z_s=bus_status(1,var1);
            bus_status(1,var1)=bus_status(1,var2);
            bus_status(1,var2)=z_s;
        end
    end
end