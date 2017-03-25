function draw_force(steps, force,lcmgl)
   for j = 1:4
    lcmgl.glPushMatrix();
    lcmgl.glTranslated(steps(j).pos(1),...
                       steps(j).pos(2),...
                       steps(j).pos(3));
    axis = [force(1:3,j)',4];
    axis = axis([4,1:3]); % LCMGL wants [angle; axis], not [axis; angle]
    axis(1) = 180 / pi * axis(1);
    lcmgl.glRotated(axis(1), axis(2), axis(3), axis(4));
    lcmgl.sphere([0;0;0], 0.015, 20, 20);
    lcmgl.glPushMatrix();
    len = 0.05;
    lcmgl.glTranslated(len / 2, 0, 0);
    lcmgl.drawArrow3d(len, 0.02, 0.02, 0.005);
    lcmgl.glPopMatrix();
    lcmgl.glPopMatrix();
    end
end