function [H] = fbmod_draw_timers(H,T)

T.xc = T.pos(:,1);
T.yc = T.pos(:,2);

for i=1:height(T)
    H(end+1).name = T.name{i};
    H(end).type = 'timer';
    switch(T.style{i})
        case {'node'}
            H(end).handles = draw_node(T(i,:));
        case {'box'}
            H(end).handles = draw_box(T(i,:));
    end
    H(end).hvec = [H(end).handles.fh H(end).handles.th];
    if ~T.has_action(i)
        set([H(end).handles.fh H(end).handles.th],'Visible','off');
    end
end

end