function [Xnew,Ynew,Znew]=AxisAngleRotate2(Xdata,Ydata,Zdata,direction,alpha)

    thetaX = alpha(1);
    thetaY = alpha(2);
    thetaZ = alpha(3);

    if direction(1)==1
        Y_x = Ydata*cos(thetaX) - Zdata*sin(thetaX);
        Z_x = Ydata*sin(thetaX) + Zdata*cos(thetaX);

        Ydata=Y_x;
        Zdata=Z_x;
    end

    if direction(2)==1
        X_y = Xdata*cos(thetaY) + Zdata*sin(thetaY);
        Z_y = Zdata*cos(thetaY) - Xdata*sin(thetaY);

        Xdata=X_y;
        Zdata=Z_y;
    end

    if direction(3)==1
        X_z = Xdata*cos(thetaZ) - Ydata*sin(thetaZ);
        Y_z = Xdata*sin(thetaZ) + Ydata*cos(thetaZ);

        Xdata=X_z;
        Ydata=Y_z;
    end

    Xnew=Xdata;
    Ynew=Ydata;
    Znew=Zdata;
end