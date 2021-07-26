function [es,fs] = steelking(Es,fy,fsu,esh,esu,dels)

r = esu - esh;
m = ((fsu/fy)*((30*r+1)^2)-60*r-1)/(15*(r^2));
es = [0:dels:esu];
ey = fy/Es;

for i=1:length(es)
    if es(i)<ey
        fs(i) = Es*es(i);
    end
    if es(i)>=ey & es(i)<=esh
        fs(i) = fy;
    end
    if es(i)>esh
        fs(i) = ((m*(es(i)-esh)+2)/(60*(es(i)-esh)+2) + (es(i)-esh)*(60-m)/(2*((30*r+1)^2)))*fy;
    end
end


