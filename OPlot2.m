function [r_e, r_p, r_h] = Oplot2(m1, m2, e, rp, th_star, w, r1 ) 
    % INPUTS: m1 [kg], m2[kg], e, rp [km], th_star [RAD], w, r1 (radius of
    % m1 for drawing purposes) [km]

    %DEFS
    w = w*ones(1,length(th_star));

    G = 6.6743*10^-20; %km^3 / s^2
    a = rp/(1-e);
    p = a*(1-e^2);

    %Part D
    %Processing

    r = (p*ones(1,length(th_star-w)))./(ones(1,length(th_star))+e*cos(th_star-w));

    r_e = r.*cos(th_star);
    r_p = r.*sin(th_star);
    r_h = zeros(1,length(th_star));
    
    circ = linspace(0,2*pi,10000);

    earthX = r1.*cos(circ);
    earthY = r1.*sin(circ);

    %Plotting
    figure()
    plot(r_e,r_p, 'c', 'DisplayName','Spacecraft Orbit');
    hold on
    %plot([0, 7150],[0, 0], 'g', 'DisplayName', 'r bar');
    hold on
    %plot([7150, 7150],[0, 8.675*scale], 'r', 'DisplayName','v bar');
    hold on 
    plot(earthX, earthY, 'r', 'DisplayName', 'Body 1');
    hold on
   

    title("Zvirani Oplot" , 'interpreter','latex');
    xlabel("e-hat direction km" , 'interpreter','latex');
    ylabel("p-hat direction km" , 'interpreter','latex'); 
    %axis([-10*10^4 4*10^4 -5*10^4 50000]) 
    axis equal
    grid on
    legend('location','northwest')
end
