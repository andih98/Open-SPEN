function im=sos(data,dir)

if nargin==2
% im=sqrt(sum(abs(data.^2),dir));
sos=abs(sum(data.*conj(data),dir)).^(1/2);
else

sos=abs(sum(data.*conj(data),3)).^(1/2);


end
im=sos./max(sos(:));  


end
