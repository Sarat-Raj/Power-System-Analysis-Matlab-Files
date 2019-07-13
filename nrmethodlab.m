                %  NR METHOD USING POLAR CO-ORDINATES FOR 14 BUS SYSTEM 
                
clc
clear all;
load bus5.m
z=bus5;
fid=fopen('output5.m','w');
n=z(1,1);            % No. of buses
m=z(1,2);	  	 	 % no. of lines
t=z(1,3);            % no. of transformers
fs=z(1,4);           % no. of fixed shunts
nv=z(1,5);           % no. of PV buses
dvlmt=z(1,6);        % convergence limit
for i=1:m+t
     Yhlc(i)=0;
     Yfsh(i)=0;
     Yrsh(i)=0;
     Yrtl(i)=0;
     Yrtr(i)=0;
     fshl(i)=0;
 end

 % OFF DIAGONAL ELEMENTS OF YBUS:
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:m+t            % FOR TRANSMISSION LINES 
     a=z(i+1,1);
     b=z(i+1,2);   
     if i<=m
         Yhlc(a)=Yhlc(a)+complex(0,z(i+1,5));
         Yhlc(b)=Yhlc(b)+complex(0,z(i+1,5));
         Ybus(a,b)=(-1)/complex(z(i+1,3),z(i+1,4));
         Ybus(b,a)=Ybus(a,b);
         Ymag(a,b)=abs(Ybus(a,b));
         Ymag(b,a)=abs(Ybus(a,b));
         Yang(a,b)=angle(Ybus(a,b));   
         Yang(b,a)=Yang(a,b);
     else      %  FOR TRANSFORMERS SB,EB,R,X,TAP
         tap=z(i+1,5);
         Ybus(a,b)=(-1)/complex(z(i+1,3),z(i+1,4));
         Yrtl(i)=(1-(1/tap))*Ybus(a,b);
         Yrsh(b)=Yrsh(b)+Yrtl(i);
         Yrtr(i)=((-1)*Yrtl(i)/tap);
         Yrsh(a)=Yrsh(a)+Yrtr(i);
         Ybus(a,b)=Ybus(a,b)/tap;
         Ybus(b,a)=Ybus(a,b);
         Ymag(a,b)=abs(Ybus(a,b));
         Ymag(b,a)=Ymag(a,b);
         Yang(a,b)=angle(Ybus(a,b));
         Yang(b,a)=Yang(a,b);    
     end
end
for i=1:fs           %  No of fixed shunts, bus number, shunt-admittance
    p=z(1+m+t+i,1);
    Yfsh(i)=1/complex(0,z(1+m+t+i,2));
    Yrsh(p)=Yrsh(p)-Yfsh(i);
end

 % DIAGONAL ELEMENTS OF YBUS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n             
     sum1=0;
     for j=1:n
         if Ybus(i,j)~=0
             sum1=sum1-Ybus(i,j);
         end
     end
     Ybus(i,i)=sum1+Yhlc(i)-Yrsh(i);
     Ymag(i,i)=abs(Ybus(i,i));
     Yang(i,i)=angle(Ybus(i,i));
end
Ybus
                      % READING OF BUS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
for i=1:n
    bn(i)=z(1+m+t+fs+i,1);
    type(i)=z(1+m+t+fs+i,2);
    typenew(i)=type(i);
    Pg(i)=z(1+m+t+fs+i,3);
    Qg(i)=z(1+m+t+fs+i,4);
    Pd(i)=z(1+m+t+fs+i,5);
    Qd(i)=z(1+m+t+fs+i,6);
    vsp(i)=z(1+m+t+fs+i,7);
    del(i)=0;
    Psch(i)=Pg(i)-Pd(i);
    Qsch(i)=Qg(i)-Qd(i);
end
for i=1:nv
    bus=z(1+m+t+fs+n+i,1);
    Vsp(bus)=z(1+m+t+fs+n+i,2);
    Qmin(bus)=z(1+m+t+fs+n+i,4);
    Qmax(bus)=z(1+m+t+fs+n+i,5);
end
for g=1:10     %---------------LOAD FLOW ITERATIONS-------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MISMATCH POWER CALUCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
if(type(i)~=3)        
    sum3=0;
    sum4=0;
    for j=1:n
            sum1=vsp(i)*vsp(j)*Ymag(i,j);
            sum2=Yang(i,j)+del(j)-del(i);
            sum3=sum3+sum1*sin(sum2);
            sum4=sum4+sum1*cos(sum2);
    end
    Pcal(i)=sum4;
    Qcal(i)=-sum3;
end
end
Qcal
for i=1:n
if typenew(i)==1
    if Qcal(i)<Qmin(i)
        Qcal(i)=Qmin(i);
        Qg(i)=Qd(i)+Qcal(i);
        Qsch(i)=Qg(i)-Qd(i);
        type(i)=0;
    elseif Qcal(i)>Qmax(i);
        Qcal(i)=Qmax(i);
        Qg(i)=Qd(i)+Qcal(i);
        Qsch(i)=Qg(i)-Qd(i);
        type(i)=0;        
    elseif Qcal(i)>=Qmin(i) & Qcal(i)<=Qmax(i)
        type(i)=1;        
    end
Qg(i)=Qd(i)+Qcal(i);
end
end
%%%%%%%%%%%%%%%%%%% DELTA P & Q 
for i=1:n
if(type(i)~=3)        
    if type(i)==1
        MM(2*i-1,1)=Psch(i)-Pcal(i);
        MM(2*i,1)=0;
    end
    if type(i)==0
        MM(2*i-1,1)=Psch(i)-Pcal(i);
        MM(2*i,1)=Qsch(i)-Qcal(i);    
    end
end
end
MM
%%%%%%%%%%%%%%%%%%% TEST FOR CONVERGENCE
limit=max(abs(MM))          
if limit<dvlmt
    break;
end
%%%%%%%%%%%%%%%%%%%%%%%% JACOBIAN  ELEMENTS %%%%%%%%%%%%%%%%%%
for i=1:n
    for j=1:n
        if (i==j)
            sum1=0;
            sum2=0;
            for k=1:n
                if (i~=k)
                    sum1=sum1+vsp(i)*vsp(k)*Ymag(i,k)*sin(Yang(i,k)+del(k)-del(i));
                    sum2=sum2+vsp(i)*vsp(k)*Ymag(i,k)*cos(Yang(i,k)+del(k)-del(i));
                end
            end
            JBN(2*i-1,2*i-1)=sum1;
            JBN(2*i-1,2*i)=sum2+2*vsp(i)^2*real(Ybus(i,i));
            JBN(2*i,2*i-1)=sum2;
            JBN(2*i,2*i)=-sum1-2*vsp(i)^2*imag(Ybus(i,i));
        else
            JBN(2*i-1,2*j-1)=(-1)*vsp(i)*vsp(j)*Ymag(i,j)*sin(Yang(i,j)+del(j)-del(i));
            JBN(2*i-1,2*j)=vsp(i)*vsp(j)*Ymag(i,j)*cos(Yang(i,j)+del(j)-del(i));
            JBN(2*i,2*j-1)=(-1)*vsp(i)*vsp(j)*Ymag(i,j)*cos(Yang(i,j)+del(j)-del(i));
            JBN(2*i,2*j)=(-1)*vsp(i)*vsp(j)*Ymag(i,j)*sin(Yang(i,j)+del(j)-del(i));
        end
    end
end

for i=1:n
    if(type(i)==3)    
        JBN(2*i-1,1:2*n)=0;
        JBN(1:2*n,2*i-1)=0;
        JBN(2*i-1,2*i-1)=1;
    end    
end
for i=1:n
    if(type(i)==1) || (type(i)==3)       
        JBN(2*i,1:2*n)=0;
        JBN(1:2*n,2*i)=0;
        JBN(2*i,2*i)=1;
    end    
end
JBN
COR=inv(JBN)*MM;

% ADDITON OF CORRECTION VECTOR ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for i=1:n
    if  (type(i)~=3)           
        del(i)=del(i)+COR(2*i-1,1);
    end
    if type(i)==0
        vsp(i)=vsp(i)*(1+COR(2*i,1));
    end
    e(i)=vsp(i)*cos(del(i))
    f(i)=vsp(i)*sin(del(i))
    Vc(i)=complex(e(i),f(i));
end

% SLACK BUS POWER CALUCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
for j=1:n
    if (type(j)==3)
        slack=j;
    end
end
Spq=0;
for j=1:n
    Spq=Spq+conj(Vc(slack))*Ybus(slack,j)*Vc(j);
end
Pcal(slack)=real(Spq);
Qcal(slack)=-imag(Spq);
Pg(slack)=Pcal(slack)+Pd(slack);
Qg(slack)=Qcal(slack)+Qd(slack);
end    %%%%%%%%%%%%%%%%%%%%%%%%  MAIN LOOP END  %%%%%%%%%%%%%%

    % OUTPUT AT THE END OF ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
fprintf(fid,'ITERATION-%i\n',g-1);
fprintf(fid,'-------------------------------------------------------------------------------------------\n');
fprintf(fid,'Bus\t\tV0ltage\t\tDelta\t\t\tPg\t\t\tQg\t\t\tPd\t\t\tQd\n');
fprintf(fid,'\t\t(p.u.)\t\t(degrees)\t\t(MW)\t\t(MVAR)\t\t(MW)\t\t(MVAR)\n');
fprintf(fid,'-------------------------------------------------------------------------------------------\n');
for i=1:n
fprintf(fid,'%i\t\t%4.3f\t\t%4.3f\t\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\n',bn(i),vsp(i),(180/pi)*del(i),Pg(i)*100,Qg(i)*100,Pd(i)*100,Qd(i)*100);
end
fprintf(fid,'-------------------------------------------------------------------------------------------\n');


% CALCULATION OF LINE FLOWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    a=vsp(i)*cos(del(i));
    b=vsp(i)*sin(del(i));
    Vc(i)=complex(a,b);
    fshl(i)=0;
end
sum3=0;
slack=0;
fprintf(fid,'\nTRANSMISSION LINE FLOWS\n');
fprintf(fid,'-------------------------------------------------------------------------------------\n');
fprintf(fid,'\t\t\tForward Power Flow\t\tReverse Power Flow\t\tPower Losses\n');
fprintf(fid,'Sb\tEb\t\tMW\t\t\tMVAR\t\tMW\t\t\tMVAR\t\tMW\t\t\tMVAR\n');
fprintf(fid,'-------------------------------------------------------------------------------------\n');
for i=1:m            %LINE FLOWS FOR TRANSMISSION LINES 
    a=z(i+1,1);
    b=z(i+1,2);
    Yhlc=complex(0,z(i+1,5));
    Y=(-1)*Ybus(a,b);
    Sij=Vc(a)*(((conj(Vc(a)))-(conj(Vc(b))))*conj(Y)+conj(Vc(a))*conj(Yhlc));
    Pij=real(Sij);
    Qij=imag(Sij);
    Sji=Vc(b)*(((conj(Vc(b)))-(conj(Vc(a))))*conj(Y)+conj(Vc(b))*conj(Yhlc));
    if a==1  |  b==1
        slack=slack+Sji;
    end
    Pji=real(Sji);
    Qji=imag(Sji);
    PLij=Sij+Sji;
    sum3=sum3+PLij;
    Ploss=real(PLij);
    Qloss=imag(PLij);
    fprintf(fid,'%i\t%i\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\n',a,b,Pij*100,Qij*100,Pji*100,Qji*100,Ploss*100,Qloss*100);   
end  
fprintf(fid,'-------------------------------------------------------------------------------------\n');
fprintf(fid,'TRANFORMERS LINE FLOWS\n');
fprintf(fid,'-------------------------------------------------------------------------------------\n');
fprintf(fid,'\t\t\tForward Power Flow\t\tReverse Power Flow\t\tPower Losses\n');
fprintf(fid,'Sb\tEb\t\tMW\t\t\tMVAR\t\tMW\t\t\tMVAR\t\tMW\t\t\tMVAR\n');
fprintf(fid,'-------------------------------------------------------------------------------------\n');
sum4=0;
for  i=m+1:m+t      %  FOR TRANSFORMERS SB,EB,R,X,TAP
    a=z(i+1,1);
    b=z(i+1,2); 
    Ytl=-Yrtr(i);
    Ytr=-Yrtl(i);
    Y=(-1)*Ybus(a,b);
    Sij=Vc(a)*(((conj(Vc(a)))-(conj(Vc(b))))*conj(Y)+conj(Vc(a))*conj(Ytl));
    Pij=real(Sij);
    Qij=imag(Sij);
    Sji=Vc(b)*(((conj(Vc(b)))-(conj(Vc(a))))*conj(Y)+conj(Vc(b))*conj(Ytr));
    if a==1  |  b==1
        slack=slack+Sji;
    end
    Pji=real(Sji);
    Qji=imag(Sji);
    PLij=Sij+Sji;
    sum4=sum4+PLij;
    Ploss=real(PLij);
    Qloss=imag(PLij);
    fprintf(fid,'%i\t%i\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\n',a,b,Pij*100,Qij*100,Pji*100,Qji*100,Ploss*100,Qloss*100);   
end
fprintf(fid,'-------------------------------------------------------------------------------------\n');
sum5=0;
for i=1:fs     %%%%%%%% SHUNT  (INDUCTIVE)  %%%%%%%%%
   p=z(1+m+t+i,1);
   FLoss(i)=Vc(p)*conj(Vc(p))*Yfsh(i);
   sum5=sum5-(FLoss(i));
end
sum6=sum3+sum4+sum5;

% OUTPUT FOR LINE LOSSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'-----------------------------------------------------------\n');
fprintf(fid,'***************   SYSTEM-GRID  TOTALS    ******************\n');
fprintf(fid,'-----------------------------------------------------------\n\n');
fprintf(fid,'Total Generation    : %4.4f MW     %4.4f MVAR\n',sum(Pg)*100,sum(Qg)*100); 
fprintf(fid,'Total T/L Loss      : %4.4f MW     %4.4f MVAR\n',real(sum3)*100,imag(sum3)*100); 
fprintf(fid,'Total T/F Loss      : %4.4f MW     %4.4f MVAR\n',real(sum4)*100,imag(sum4)*100); 
fprintf(fid,'Shunt (inducive)    :              %4.4f MVAR\n',imag(sum5)*100); 
fprintf(fid,'Total P - Q  Load   : %4.4f MW     %4.4f MVAR\n',sum(Pd)*100,sum(Qd)*100);
fprintf(fid,'Total Power Losses  : %4.4f MW     %4.4f MVAR\n',real(sum6)*100,imag(sum6)*100); 
fprintf(fid,'-----------------------------------------------------------\n');
fclose(fid);