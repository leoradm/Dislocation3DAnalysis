function booleanOut = withinRange(Boundary,bNum,Point)
    % define the range for the point
    xR = Boundary(bNum).xRange;
    yR = Boundary(bNum).yRange;
    zR = Boundary(bNum).zRange;
    
    % is the point contained in this range?
    booleanOut = 0;
    if Point(1) >= xR(1) && Point(1) <= xR(2)
        if Point(2) >= yR(1) && Point(2) <= yR(2)
            if Point(3) >= zR(1) && Point(3) <= zR(2)
                booleanOut = 1;
            end
        end
    end
end
