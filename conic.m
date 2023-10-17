
classdef conic

    

    %class for creating an object and relevant methods for conic orbits
    %conic(a,e,w,O,i, main_body)

    properties
        %inputs
        a %semimajor axis (km)
        e %eccentricity
        w % argument of perigee (RAD input, rad stored)
        O % RAAN/longitude of asending node (RAD input, rad stored)
        i % inclination (RAD input, rad stored)
        r1 %radius of central body of orbit

        %non-inputs
        rp ;%periapsis
        ra ;%apoapsis [km]
        p ;% semi-latus rectum [km]
        b ;% semi-minor axis
        T; %orbital period
        n ; %mean motion
        mu;
     

    end

    methods
        function obj = conic(a,e,w,O,i, main_body)
            %bodyDict
            rDict = dictionary("earth", 6378.14);
            muDict = dictionary("earth", 398600.4415);
            obj.a = a;
            obj.e = e;
            obj.w = w;
            obj.O = O;
            obj.i = i;
            obj.r1 = rDict(main_body);
            obj.mu = muDict(main_body);
            obj.n = sqrt(obj.mu/obj.a^3);

            if obj.e>1
                obj.rp = abs(obj.a)*(obj.e-1);
                obj.ra = abs(obj.a)*(obj.e+1);
                obj.p = abs(obj.a)*((obj.e^2)-1);
                obj.b = abs(obj.a)*sqrt((obj.e^2)-1);
                obj.T = 0;


            elseif obj.e<1
                obj.rp = obj.a*(1-obj.e);
                obj.ra = obj.a*(1+obj.e);
                obj.p = obj.a*(1-obj.e^2);
                obj.b = obj.a*sqrt(1-obj.e^2);
                obj.T = 2*pi*sqrt((obj.a^3/obj.mu));
            else
            end

        end

        function [r_e, r_p, r_h] = plot2(obj, th1, th2)

            th_star = linspace(deg2rad(th1), deg2rad(th2), 1000);

            %r = (obj.p*ones(1,length(th_star-obj.w)))./(ones(1,length(th_star))+obj.e*cos(th_star-obj.w));
            r = (obj.p*ones(1,length(th_star)))./(ones(1,length(th_star))+obj.e*cos(th_star));
            r_e = r.*cos(th_star);
            r_p = r.*sin(th_star);
            r_h = zeros(1,length(th_star));

            circ = linspace(0,2*pi,100000);

            earthX = obj.r1.*cos(circ);
            earthY = obj.r1.*sin(circ);


            %Plotting
            figure()
            plot(r_e,r_p, 'r', 'DisplayName','Spacecraft Orbit');
            hold on
            %plot([0, 7150],[0, 0], 'g', 'DisplayName', 'r bar');
            hold on
            %plot([7150, 7150],[0, 8.675*scale], 'r', 'DisplayName','v bar');
            hold on
            plot(earthX, earthY, 'b', 'DisplayName', 'Body 1');
            hold on


            title("Zvirani 2D-Plot, Orbit Normal" , 'interpreter','latex');
            xlabel("e-hat direction km" , 'interpreter','latex');
            ylabel("p-hat direction km" , 'interpreter','latex');
            %axis([-10*10^4 4*10^4 -5*10^4 50000])
            axis equal
            grid on
            legend('location','northwest')
        end

        function sc_pos = plot3(obj, th0, th1, suppress)

            % PURPOSE: Function to plot a given conic section, two-body orbit in 3D.
            % INPUTS: th0, th1: initial and final true anomaly through which to plot [deg],  suppress: if 1, plots will be suppressed.
            %OUTPUT: length(t) x 3 (X,Y,Z) array of spacecraft position in inertial equatorial Xhat,
            %Yhat, Zhat coordinates

            %input cleaning
            th_star = linspace(deg2rad(th0),deg2rad(th1),1000);




            %DEFS
            local_w = obj.w*ones(1,length(th_star));
            RAAN = obj.O*ones(1,length(th_star));
            local_i = obj.i*ones(1,length(th_star));
            th = th_star - local_w;

            G = 6.6743*10^-20; %km^3 / s^2


            X_vec_vis = 5*10^4;
            Y_vec_vis = 5*10^4;
            Z_vec_vis = 5*10^4;



            %Processing

            r = (obj.p*ones(1,length(th)))./(ones(1,length(th))+obj.e*cos(th));
            h_coord = zeros(1,length(r));


            dcm = angle2dcm(RAAN, local_i, th_star, 'ZXZ');

            sc_pos = zeros(length(th),3);

            for local_i =1:length(th_star)
                sc_pos(local_i,:) = dcm(:,:,local_i)'*[r(local_i),0,0]';

            end


            for local_i = 1:length(sc_pos)
                mag = sqrt(sc_pos(local_i,1)^2 + sc_pos(local_i,2)^2 + sc_pos(local_i,3)^2);
                if mag<= obj.r1
                    fprintf("Spacecraft Impact Detected")
                end
            end


            [Ex,Ey,Ez] = sphere(50);
            if suppress ~= 1
                figure()

                plot3(sc_pos(:,1),sc_pos(:,2),sc_pos(:,3), 'r')
                hold on
                quiver3(0,0,0, X_vec_vis,0,0,'g')
                hold on
                quiver3(0,0,0, 0, Y_vec_vis,0,'g')
                hold on
                quiver3(0,0,0, 0,0, Z_vec_vis,'g')
                hold on
                set(gcf,'color','w')
                surf(Ex*obj.r1, Ey*obj.r1, Ez*obj.r1)
                colormap summer
                shading interp
                hold on
                grid on
                axis equal
                title("Spacecraft Orbit", 'interpreter','latex')
                xlabel('X (km)', 'interpreter','latex')
                ylabel('Y (km)', 'interpreter','latex')
                zlabel('Z (km)', 'interpreter','latex')
                %legend( 'location','northwest')

                %         xlim([-X_vec_vis X_vec_vis])
                %         ylim([-Y_vec_vis*3 Y_vec_vis])
                %         zlim([-Z_vec_vis Z_vec_vis])
                view(3)


            end




        end


    end




end



