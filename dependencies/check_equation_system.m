function [] = check_equation_system(min_number, is_number, grid_type)
% grid_type: 'volume' or 'surface'

if strcmp(grid_type, 'volume')
    string = '';
elseif strcmp(grid_type, 'surface')
    string = ' pairs of';
end
    
% make sure that we're getting an overdetermined equation system
fprintf('The minimum required number of%s sampling points for the desired N is %d. %d are available.\n\n', string, min_number, is_number);
    
% if underdetermined
if is_number < min_number
    warning('The equation system is underdetermined (%d required vs. %d available).', min_number, is_number);
end
    

end

