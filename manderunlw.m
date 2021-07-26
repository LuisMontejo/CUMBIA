
function [ec,fcu] = manderunlw(Ec,nbl,Dbl,Dh,clb,s,fpc,fyh,eco,esm,espall,section,D,d,b,ncx,ncy,wi,dels)


% unconfined concrete:

ec = [0:dels:espall];
Esecu = fpc/eco;
ru    = Ec/(Ec-Esecu);
xu    = ec./eco;
ru2   = Ec/(Ec-1.8*fpc/eco);

for i = 1:length(ec)
    if ec(i)<eco
        fcu(i) = fpc*xu(i)*ru/(ru-1+xu(i)^ru);
    end
    if ec(i)>=eco & ec(i)<1.3*eco
        fcu(i) = fpc*xu(i)*ru2/(ru2-1+xu(i)^ru2);
    end
    if (ec(i)>=1.3*eco & ec(i)<=espall)
        fcu(i) = fpc*(1.3*ru2/(ru2-1+1.3^ru2))*(1-(ec(i)-1.3*eco)/(espall-1.3*eco));
    end
    if ec(i)>espall
        fcu(i) = 0;
    end
end

