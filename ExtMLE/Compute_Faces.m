function Faces = Compute_Faces(Complex,Nvar)
% Compute all the faces of the Simplicial Complex
% Uses the function sub_by_size_next to get all subset of a set.
n=length(Complex);
Faces{1}=0;
for i=1:Nvar
    Faces{i+1}=i;
end

for i=1:n % Loop trhought Complex 
    a_new=[];
    size_new=0;
    more_new=0;
    lenFace = length(Complex{i});
    while((size_new>1)||(size_new==0))
        % First call
        [a_new, size_new, more_new] = sub_by_size_next(lenFace,a_new,size_new,more_new);
        % Check whether to add to Faces or not
        if(size_new>1)
            add=1;
            for j=(Nvar+2):length(Faces)
                if(isequal(Complex{i}(a_new), Faces{j})) % If face already present, break
                    add=0;
                    break;
                end             
            end
            if(add==1) % Face is new: add to the list of faces
               Faces{length(Faces)+1}=Complex{i}(a_new);
            end             
        end
    end
end

