function Int=Make_Integer(Y)
% Make_Integer
%   Convert rational matrix Y into an integer matrix.
% Inputs:
%   Y: rational matrix.
% Outputs:
%   Int: integer matrix Int so that Y = Int ./ lstcommult.
[num,den]=rat(Y);
dens=sort(unique(abs(den)));
lstcommult=1;
for i=1:length(dens)
    if(dens(i)>1)
       lstcommult=lcm(lstcommult,dens(i));
    end
end
%Int = num./den.*lstcommult;
Int = Y.*lstcommult;

