function [] = display_progress(in)
% in: either label string or relative progress

persistent label;

display_width = 50;

if isempty(label) && ischar(in)
    
    fprintf('%s\n', in);
    label = -1;

elseif ~isempty(label) && ischar(in)

    label = [];  
    fprintf([in '\n']);

elseif isnumeric(in)
        
    nDots = floor(in * display_width);
    str_out = ['|' repmat('-', 1, nDots-1) '*' repmat('-',1,display_width-nDots) '|'];

    % print progress
    if label == -1
        % don't do carriage return during first run
        fprintf(str_out);
    else
        % do it during all the other runs
        fprintf([label str_out]);
    end
    
    % Update carriage return
    label = repmat('\b',1,length(str_out));
    
end

