function [es,fs] = Raynor(Es,fy,fsu,esh,esu,dels,C1,Ey)


es  = [0:dels:esu];
ey  = fy/Es;
fsh = fy + (esh-ey)*Ey

for i=1:length(es)
    if es(i)<ey
        fs(i) = Es*es(i);
    end
    if es(i)>=ey & es(i)<=esh
        fs(i) = fy+(es(i)-ey)*Ey;
    end
    if es(i)>esh
        fs(i) = fsu-(fsu-fsh)*(((esu-es(i))/(esu-esh))^C1);
    end
end


