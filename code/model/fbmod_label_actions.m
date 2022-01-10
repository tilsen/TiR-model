function [h] = fbmod_label_actions(T,M)

h=[];
for i=1:height(T)
   
    targ_ixs = T.trg{i};
    for j=1:length(targ_ixs)
        
        force = T.F{i}(:,j);
        valence = T.chi{i}(j);
        targ_name = M.S.name{targ_ixs(j)};
        
        if contains(targ_name,'gate') && ~contains(T.label{i},'theta')
            gest_id = M.S.id(ismember(M.S.name,strrep(targ_name,'gate_','gest_')));   
            targ_lab = M.G.label{ismember(M.G.id,gest_id)};
            str = [T.label{i} '{\rightarrow}' targ_lab];
        else
            %!!not implemented
            str = [T.label{i}];
        end
        
        str = regexprep(str,'(?<!^)\$(?!$)','');
        
        t_action = M.t(find(force~=0,1,'first'));
                
        switch(valence)
            case -1
                va = 'top';
            case 1
                va = 'bot';
        end
        
        sfac = 1.1;
        
        if ~isempty(t_action)
            h(end+1) = text(t_action,valence*sfac,str,...
                'verti',va,'hori','center','interp','latex');
        end
                
        
    end
    
end



end