function blue_white_red = blueWhiteRedColormap(Nc)


% Nc = 256;
blue_white_r = linspace(0,1,Nc/2);
blue_white_g = linspace(0,1,Nc/2);
blue_white_b = ones(1,Nc/2);
blue_white = [blue_white_r', blue_white_g', blue_white_b'];

white_red_r = ones(1,Nc/2);
white_red_g = linspace(1,0,Nc/2);
white_red_b = linspace(1,0,Nc/2);
white_red = [white_red_r', white_red_g', white_red_b'];
blue_white_red = [blue_white; white_red];



