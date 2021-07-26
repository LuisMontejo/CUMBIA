
function [ec,fcu] = manderun(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,section,D,d,b,ncx,ncy,wi,dels)

% unconfined concrete:

ec = [0:dels:espall];
Esecu = fpc/eco;
ru    = Ec/(Ec-Esecu);
xu    = ec./eco;

for i = 1:length(ec)
    if ec(i)<2*eco
        fcu(i) = fpc*xu(i)*ru/(ru-1+xu(i)^ru);
    end
    if (ec(i)>=2*eco & ec(i)<=espall)
        fcu(i) = fpc*(2*ru/(ru-1+2^ru))*(1-(ec(i)-2*eco)/(espall-2*eco));
    end
    if ec(i)>espall
        fcu(i) = 0;
    end
end

