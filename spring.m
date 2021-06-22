function [xs ys] = spring(xa,ya,xb,yb,varargin)
persistent ne Li_2 ei b
if nargin > 4 
    [ne a r0] = varargin{1:3};                  
    Li_2 =  (a/(4*ne))^2 + r0^2;                
    ei = 0:(2*ne+1);                            
    j = 0:2*ne-1; b = [0 (-ones(1,2*ne)).^j 0]; 
end
R = [xb yb] - [xa ya]; mod_R = norm(R); 
L_2 = (mod_R/(4*ne))^2; 
if L_2 > Li_2
   error('Spring:TooEnlargement', ...
   'Initial conditions cause pulling the spring beyond its maximum large. \n Try reducing these conditions.')
else
    r = sqrt(Li_2 - L_2);   
end
c = r*b;    
u1 = R/mod_R; u2 = [-u1(2) u1(1)];  
xs = xa + u1(1)*(mod_R/(2*ne+1)).*ei + u2(1)*c; 
ys = ya + u1(2)*(mod_R/(2*ne+1)).*ei + u2(2)*c; 