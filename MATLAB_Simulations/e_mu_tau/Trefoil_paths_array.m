% Trefoil_paths_array

clear
clc

N_loops = 1;
t = 0:0.02:N_loops*2*pi;
R = 10;
array_max1 = 4;
array_max2 = 4;

make_test = 1;
if make_test
    
    N_loops = 1;
    t = 0:0.02:N_loops*2*pi;
    nr1 = 2;
    nr2 = nr1/1.5;
    nz = 3;
    x = R*(nr2*cos(1*t) - nr1*cos(nr1*t));
    y = R*(nr2*sin(1*t) + nr1*sin(nr1*t));
    z = R*sin(nz*t) + 0*t;
    % Original:
%     x = R*(cos(t) - nr*cos(nr*t));
%     y = R*(sin(t) + nr*sin(nr*t));
%     z = -nz*R*sin(nz*t);

%     x = R*(2 + cos(10*t)).*cos(2*t);
%     y = R*(2 + cos(10*t)).*sin(2*t);
%     z = R*sin(3*t);

%     x = R*(2 + cos(3*t)).*cos(2*t);
%     y = R*(2 + cos(3*t)).*sin(2*t);
%     z = R*sin(3*t);

    use_coil_on_coil = 1;
    if use_coil_on_coil
        [theta,phi] = meshgrid(linspace(0,2*pi,50));
        r = 2;
        R = r*5;
        n = 10;
        v = 0;
        
%         x = (R + r*cos(theta)).*cos(phi);
%         y = (R + r*cos(theta)).*sin(phi);
%         z = r*sin(theta);
        
        xs = (R + r*cos(theta)).*cos(phi);
        ys = (R + r*cos(theta)).*sin(phi);
        zs = r*sin(theta);
        
        X = R*sin(t);
        Y = R*cos(t);
        Z = v*t;
        
        x = X + r.*cos(n*t).*cos(t);
        y = Y + r.*cos(n*t).*sin(t);
        z = Z + r.*sin(n*t);
        
        figure(9)
        c_alpha = 0.5;
%         plot3(X,Y,Z,'-')
        surf(xs,ys,zs,'FaceAlpha',c_alpha)
        shading interp
        colormap gray
        
        hold on
        plot3(x,y,z,'-','LineWidth',2)
        hold off
        
        grid on
        axis equal
    end
    
    figure(10)
    subplot(1,2,1)
    plot3(x,y,z,'-')
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['nr = ',num2str(nr1),', nz = ',num2str(nz)])
    
    subplot(1,2,2)
    plot3(x,y,z,'-')
    view(2)
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['nr = ',num2str(nr1),', nz = ',num2str(nz)])
end
%
%% 3D view:
if ~make_test
    figure(11)
    for ind1 = 1:array_max1
        for ind2 = 1:array_max2
            nr1 = ind1;
            nz = ind2;

            x = R*(cos(t) + nr1*cos(nr2*t));
            y = R*(sin(t) + nr1*sin(nr2*t));
            z = R*sin(nz*t);


            subplot(array_max1,array_max1,(ind2-1)*(array_max2) + ind1)
            plot3(x,y,z,'-')
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title(['nr = ',num2str(nr1),', nz = ',num2str(nz)])

        end
    end
%
%% 2D view:
    figure(12)
    for ind1 = 1:array_max1
        for ind2 = 1:array_max2
            nr1 = ind1;
            nz = ind2;

            x = R*(cos(t) + nr1*cos(nr1*t));
            y = R*(sin(t) + nr1*sin(nr1*t));
            z = R*sin(nz*t);


            subplot(array_max1,array_max1,(ind2-1)*(array_max2) + ind1)
            plot3(x,y,z,'-')
            view(2)
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title(['nr = ',num2str(nr1),', nz = ',num2str(nz)])

        end
    end
end
%
