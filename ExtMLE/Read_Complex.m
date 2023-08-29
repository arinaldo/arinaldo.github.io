function Table = Read_Complex(filename)
% Read in the table and model quantities
% Input: filename
% OUtput: Structure

fid=fopen(filename,'r');
% Number of variables
Nvar=fscanf(fid,'%g',1);
% Vector of length Nvar with number of categories
Ncat=fscanf(fid,'%g',Nvar);% Row vector of length Nvar
% Hierarchical or not
%Hier=textscan(fid,'%c',1);
Hier=fscanf(fid,'%s',1);
% Number of facet of the simplicial complex
Nfacet=fscanf(fid,'%g',1); % same as length(Table.Complex)
% Allocate memory for the complex
Complex=cell([Nfacet 1]);
% Read in complex values
temp=fscanf(fid,'%g',1);
for i=1:Nfacet
    Complex{i,1}=fscanf(fid,'%g',temp);
    temp=fscanf(fid,'%g',1);
end
% close the file
fclose(fid);

% Create the structure
Table.Nvar=Nvar;
Table.Ncat=Ncat;
Table.Hier=Hier;
Table.Complex=Complex;
    
%temp=textscan(fid,'%n',1);
%Nvar=temp{1,1};

%temp=textscan(fid,'%n',Nvar);
%Ncat=temp{1,1}; 

%temp=textscan(fid,'%s',1);
%Hierarchical=temp{1,1};
%temp=textscan(fid,'%n',1);
%Nfacet=temp{1,1};
%Complex{i,1}=textscan(fid,'%n',temp);
%Complex{i,1}=Complex{1,1}
