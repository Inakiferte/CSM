function [J,N,dN_dxhi] = shape(e,xhi,n_G,x1,x2,x3)


if n_G == 2
   
   dN_dxhi    = zeros(1,n_G);
   N          = zeros(1,n_G);
   dN_dxhi(1) = -0.5;
   dN_dxhi(2) = 0.5;
   N(1)       = 0.5 * (1 - xhi);
   N(2)       = 0.5 * (1 + xhi);

   J          = dN_dxhi(1) * x1 + dN_dxhi(2) * x2; 
end

if n_G == 3
   
   dN_dxhi    = zeros(1,n_G);
   N          = zeros(1,n_G);

   dN_dxhi(1) = xhi - 0.5;
   dN_dxhi(2) = -2 * xhi;
   dN_dxhi(3) = xhi + 0.5;
   N(1)       = 0.5 * xhi * (xhi - 1);
   N(2)       = (1 - xhi^2);
   N(3)       = 0.5 * xhi * (xhi + 1);
  
   J          = dN_dxhi(1) * x1 + dN_dxhi(2) * x2 + dN_dxhi(3) * x3;  
end

end