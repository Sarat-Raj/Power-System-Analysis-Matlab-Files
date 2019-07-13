function[Y]=zbussymfault(zdata1)
nl=zdata1(:,1);
nr=zdata1(:,2);
r=zdata1(:,3);
x=zdata1(:,4);
nbr=length(zdata1(:,1));
nbus=max(max(nl),max(nr));
z=r+1i*x;
y=ones(nbr,1)./z;
Y=zeros(nbus,nbus);
for k=1:nbr;
if nl(k)>0 && nr(k)>0
Y(nl(k),nr(k))=Y(nl(k),nr(k))-y(k);
Y(nr(k),nl(k))=Y(nl(k),nr(k));
end
end
for n=1:nbus
for k=1:nbr
if (nl(k))== n || (nr(k))== n
Y(n,n)=Y(n,n)+y(k);
else end
end
end
kv=input('Enter base KV : ');
mva=input('Enter base MVA : ');
bc=((mva*1000)/(sqrt(3)*kv));
Zbus=inv(Y);
disp('BUS IMPEDANCE MATRIX : ');
disp(Zbus);
fbn=input('Enter the Fault bus no ---> ');