
function [ec,fc] = manderconf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,espall,section,D,d,b,ncx,ncy,wi,dels,type)


% confined concrete:

sp  = s - Dh;
Ash = 0.25*pi*(Dh^2);


switch lower(section)
    case 'rectangular'
        bc   = b - 2*clb + Dh;
        dc   = d - 2*clb + Dh;
        Asx  = ncx*Ash;
        Asy  = ncy*Ash;
        Ac   = bc*dc;
        rocc = Ast/Ac;
        rox  = Asx/(s*dc);
        roy  = Asy/(s*bc);
        ros  = rox + roy;
        ke   = ((1 - sum(wi.^2)/(6*bc*dc)) * (1 - sp/(2*bc)) ...
            * (1-sp/(2*dc))) / (1 - rocc);
        ro   = 0.5*ros;
        fpl  = ke*ro*fy;

    case 'circular'
        ds   = D - 2*clb + Dh;
        ros  = 4*Ash/(ds*s);
        Ac   = 0.25*pi*(ds^2);
        rocc = Ast/Ac;
        switch lower(type)
            case 'spirals'
                ke   = (1-sp/(2*ds))/(1-rocc);
            case 'hoops'
                ke   = ((1-sp/(2*ds))/(1-rocc))^2;
            otherwise
                disp('tranverse reinforcement should be spirals or hoops'); return;
        end
        fpl  = 0.5*ke*ros*fy;

    otherwise
        disp('section not available'); return;
end


fpcc = (-1.254 + 2.254*sqrt(1 + 7.94*fpl/fpc) - 2*fpl/fpc)*fpc;
ecc  = eco*(1 + 5*(fpcc/fpc-1));
Esec = fpcc/ecc;
r    = Ec/(Ec-Esec);
ecu  = 1.5*(0.004 + 1.4*ros*fy*esm/fpcc);


ec = [0:dels:ecu];
x  = (1/ecc)*ec;
fc = fpcc*x*r./(r-1+x.^r);




