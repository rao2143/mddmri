function [Cx,Cy,Cz,div] = curldiv(Vx,Vy,Vz,dx,dy,dz)[Vxdx,Vxdy,Vxdz] = gradient(Vx,dx,dy,dz);[Vydx,Vydy,Vydz] = gradient(Vy,dx,dy,dz);[Vzdx,Vzdy,Vzdz] = gradient(Vz,dx,dy,dz);Cx = Vydz - Vzdy;Cy = Vzdx - Vxdz;Cz = Vxdy - Vydx;div = Vxdx + Vydy + Vzdz;