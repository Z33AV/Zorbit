function sc_pos = Oplot3(m1, r1, m2, e, rp, th_star, w, i, RAAN, suppress)


    % PURPOSE: Function to plot a given conic section, two-body orbit in 3D.
    % INPUTS: M1: mass of body 1, r1: radius of body 1 (for vis and impact
    % detection) [km], m2: mass of body 2, e: eccentricity of orbit, rp:
    % periapsis of orbit [km], th_star: list of true anomaly values to plot
    % through [RAD], w: argument of perigee[deg], i: inclination [deg], RAAN: Right angle of
    % ascending node [deg], suppress: if 1, plots will be suppressed.
    %OUTPUT: length(t) x 3 (X,Y,Z) array of spacecraft position in inertial equatorial Xhat,
    %Yhat, Zhat coordinates
    
    %input cleaning
    RAAN = deg2rad(RAAN);
    i = deg2rad(i);
    %th_star = linspace(0,2*pi,1000);
    e = e;
    rp = rp;
    r1 = r1;
    w = deg2rad(w);
    
    
    
    
    %DEFS
    w = w*ones(1,length(th_star));
    RAAN = RAAN*ones(1,length(th_star));
    i = i*ones(1,length(th_star));
    th = th_star - w;
    
    G = 6.6743*10^-20; %km^3 / s^2
    a = rp/(1-e);
    p = a*(1-e^2);
    
    
    X_vec_vis = 5*10^4;
    Y_vec_vis = 5*10^4;
    Z_vec_vis = 5*10^4;
    
    
    
    %Processing
    
    r = (p*ones(1,length(th)))./(ones(1,length(th))+e*cos(th));
    h_coord = zeros(1,length(r));
    
    
    dcm = angle2dcm(RAAN, i, th_star, 'ZXZ'); 
    
    sc_pos = zeros(length(th),3);
    
    for i =1:length(th_star)
       sc_pos(i,:) = dcm(:,:,i)'*[r(i),0,0]';
    
    end
    
    
    for i = 1:length(sc_pos)
        mag = sqrt(sc_pos(i,1)^2 + sc_pos(i,2)^2 + sc_pos(i,3)^2);
        if mag<= r1
            fprintf("Spacecraft Impact Detected")
        end
    end
    
    
    [Ex,Ey,Ez] = sphere(50);
    if suppress ~= 1
        figure(1)
        
        plot3(sc_pos(:,1),sc_pos(:,2),sc_pos(:,3), 'r')
        hold on
        quiver3(0,0,0, X_vec_vis,0,0,'g')
        hold on
        quiver3(0,0,0, 0, Y_vec_vis,0,'g')
        hold on
        quiver3(0,0,0, 0,0, Z_vec_vis,'g')
        hold on
        set(gcf,'color','w')
        surf(Ex*r1, Ey*r1, Ez*r1)
        hold on
        grid on
        axis equal
        title("Spacecraft Orbit", 'interpreter','latex')
        xlabel('X (km)', 'interpreter','latex')
        ylabel('Y (km)', 'interpreter','latex')
        zlabel('Z (km)', 'interpreter','latex')
        legend( 'location','northwest')
        
%         xlim([-X_vec_vis X_vec_vis])
%         ylim([-Y_vec_vis*3 Y_vec_vis])
%         zlim([-Z_vec_vis Z_vec_vis])
        view(3)

        
    end
end
