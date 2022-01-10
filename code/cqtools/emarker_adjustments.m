function [] = emarker_adjustments(H)

arrayfun(@(c)correct_point_ypos(H.emark(c),H.emark_levelypos(c),H.emark(c).Parent),1:length(H.emark));
arrayfun(@(c)correct_text_ypos(H.etext(c),H.emark(c),H.etext(c).Parent),1:length(H.etext));

end

