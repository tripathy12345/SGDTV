function yms=EWT_Meyer_Scaling(w1,gamma,N)

Mi=floor(N/2);
w=fftshift((0:2*pi/N:2*pi-2*pi/N))';
w(1:Mi)=-2*pi+w(1:Mi);

aw=abs(w);
yms=zeros(N,1);

an=1/(2*gamma*w1);
pbn=(1+gamma)*w1;
mbn=(1-gamma)*w1;

for k=1:N
   if (aw(k)<=mbn)
       yms(k)=1;
   elseif ((aw(k)>=mbn) && (aw(k)<=pbn))
       yms(k)=cos(pi*EWT_beta(an*(aw(k)-mbn))/2);
   end
end
yms=ifftshift(yms);