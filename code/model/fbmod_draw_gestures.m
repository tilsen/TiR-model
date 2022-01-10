function [H] = fbmod_draw_gestures(H,G)

G.xc = G.pos(:,1);
G.yc = G.pos(:,2);

for i=1:height(G)
   H(end+1).name = G.name{i};
   H(end).type = 'gesture';
   H(end).handles = draw_node(G(i,:));
   H(end).hvec = [H(end).handles.fh H(end).handles.th];
end



end