function [H] = fbmod_draw_sens(H,Y)

Y.xc = Y.pos(:,1);
Y.yc = Y.pos(:,2);

for i=1:height(Y)   
   H(end+1).name = Y.name{i};
   H(end).type = 'sens';
   H(end).handles = draw_box(Y(i,:),'linestyle',':');
   H(end).hvec = [H(end).handles.fh H(end).handles.th];
   if ~Y.has_action(i)
       set(H(end).hvec,'Visible','off'); 
   end
end

end