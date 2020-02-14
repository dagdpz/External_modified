function [varagout] = process_varargin(the_keys,default_values,new_varargin)

tmp = [];
n_key = length(the_keys);
n_def = length(default_values);
nvar  = length(new_varargin);

%  make sure we have the correct number
if (n_def ~= n_key)
    fprintf(2,'process_varargin: # of keys ~= # of default values %3d %3d\n',n_key,n_def);
    error('default values error');
end

%  set up the default values
for i=1:n_key
    tmp.(the_keys{i}) = default_values{i};
end

%  user overrides ?
i = 1;
j = 2;
while (i < nvar)
    if  any(ismember(the_keys,new_varargin{i})); %% !!!!!
        tmp.(new_varargin{i}) = new_varargin{j};
    end
    i = i + 2;
    j = j + 2;
end
varagout = tmp;
