function [ncell_array]=nest_cell_array(cell_array, type)

%This function accepts as input an m x n cell array. It returns an 1 x m
%cell array of 1 x n cells or doubles, as specified by the type input argument 
if (nargin < 2)
    type= 'cell';
end


if  strcmp(type , 'cell')
    ncell_array=cell(1,size(cell_array,1));
    for i=1:size(cell_array,1)
        tmpcell=cell(1,size(cell_array,2));
        for ii=1:size(cell_array,2)
            tmpcell{ii}=cell_array{i,ii};
        end
        ncell_array{1,i}=tmpcell;
        
    end

elseif  strcmp(type ,'double')
    ncell_array = cell(1,size(cell_array,1));
    for i = 1:size(cell_array,1)
        tmpdouble = zeros(1,size(cell_array,2));
        for ii = 1:size(cell_array,2)
            tmpdouble(ii) = cell_array{i,ii};
        end
        ncell_array{1,i} = tmpdouble;
    end
    
else
    disp('Error: invalid type specified')
        
end
     
    