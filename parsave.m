function [] = parsave(sxy,x,y,k,kappas,kappat,stress,sigmas,sigmat,sv,vt,vn,sphi,phi,sGphi,GPhi)
fx = flip(x);
fy = flip(y);
fks = flip(kappas);
fkt = flip(kappat);
fss = flip(sigmas);
fst = flip(sigmat);
fvt = flip(vt);
fvn = flip(vn);
fphi = flip(phi);
fGPhi = flip(GPhi);
savefile(sxy,[[-x;fx],[y;fy]]);
savefile(k,[[kappas;fks],[kappat;fkt]]);
savefile(stress,[[sigmas;fss],[sigmat;fst]]);
savefile(sv,[[vt;fvt],[vn;fvn]]);
savefile(sphi,[phi;fphi]);
savefile(sGphi,[GPhi;fGPhi]);
end
