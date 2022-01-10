function [ix_on,ix_off,t_on,t_off] = fbmod_gestural_activation(X,t)

for i=1:size(X,2)
    ix_on{i} = find(diff(X(:,i))==1)';
    ix_off{i} = find(diff(X(:,i))==-1)';
    
    if X(end,i)==1
        ix_off{i}(end+1) = size(X,1);
    end
end

t_on = cellfun(@(c){t(c)'},ix_on);
t_off = cellfun(@(c){t(c)'},ix_off);

end