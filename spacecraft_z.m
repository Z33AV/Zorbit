classdef spacecraft_z

    properties
        orbit; %conic object class input - input 0 if orbit elements are not known
        state; %3x2 matrix of [r,v] where r and v are coloumn vectors
        basis; %string indicating which vector basis state vectors are in. Options: 'xyz','rth', 'eph'
        E; %current Eccentric anomaly, input a string if not applicable
        H; %initial Hyperbolic anomaly, same deal as E
        TA; %initial true anomaly
        h; %current angular momentum
        gamma; % current flight path angle
        M; %current mean anomaly
        body; %string of name of main body
        mu;



        fgmat; %2x2 matrix of form [f,g;fdot,gdot]
    end

    methods

        function obj = spacecraft_z(orbit, state,basis,body)
            obj.orbit = orbit;
            obj.state = state;
            obj.basis = basis;
            obj.body = body;
            muDict = dictionary("earth", 398600.4415);
            obj.mu = muDict(obj.body);
        end


        function [a,e,w,O,i,E,h,gamma,M] = kepels(obj)
            %Implementation of Vallado RV2COE Algo9, ASSUMES XYZ basis for now
            %function to determine keplerian elements using state vectors,
            %also calculates current anomalies

            rmag = norm(obj.state(:,1));
            vmag = norm(obj.state(:,2));
            hbar = cross(obj.state(:,1),obj.state(:,2));
            hmag = norm(hbar);
            nbar = cross([0,0,1]',hbar);
            ebar = (1/obj.mu)*(vmag^2 - obj.mu/rmag)*obj.state(:,1) - (1/obj.mu)*dot(obj.state(:,1),obj.state(:,2))*obj.state(:,2);
            e = norm(ebar);
            eps = (vmag^2)/2  - obj.mu/rmag;

            if e~=1.0
                a = -obj.mu/(2*eps);
            else
                a = inf;
            
            end
            i = acos (dot(hbar,[0,0,1])/hmag);

            O = acos(dot(nbar,[1,0,0])/norm(nbar));
            if dot(nbar,[0,1,0])<0
                O = 2*pi-O;
            end

            w = acos( dot(nbar, ebar)/( norm(nbar)*norm(ebar) ) );
            if dot(ebar,[0,0,1])<0
                w = 2*pi-w;
            end



            E = acos( ( a - norm(obj.state(:,1)) )/(a*e));
            h  = hbar;
            gamma = acos( norm(obj.h)/( rmag*vmag ) );
            M = Mfcn(obj.E,e);


            %obj.orbit = conic(a,e,w,O,i,obj.body) - doesn't work because
            %fuck matlab OO

        end

        function out = fgcalc(obj, E_next, dt)
            %function to determine f and g and derivative values based on
            %given information. NOTE: requires a conic object to be defined
            %or kepels to be run first.

            if obj.orbit.e<1
                f = 1-( obj.orbit.a/norm(obj.state(:,1)) )*( 1-cos(E_next-obj.E) )
                g = dt + ( sin(E_next - obj.E) - (E_next - obj.E) )/obj.orbit.n
                rn = f*obj.state(:,1) + g*obj.state(:,2) 
                r_new =norm( rn )
                fdot = ( (-obj.orbit.n*obj.orbit.a^2)/(r_new * norm( obj.state(:,1) ) )) * sin(E_next - obj.E);
                gdot = 1 - (obj.orbit.a/r_new)*(1-cos(E_next-obj.E));

                obj.fgmat = [f,g;fdot,gdot];
                out = obj.fgmat

            end

        end

        function  nstate = impulse(obj,dV, alpha)
            %inputs: deltaV magnitude [km/s], Alpha angle defined as angle between
            %initial velocity (current spacecraft state vector 2) and
            %desired velocity in radians - currently BROKEN - doesnt take
            %into account zhat velo change and breaks for small values of O
            %and w
            
            v_minus_dir = (1/norm(obj.state(:,2)))*obj.state(:,2);
            dV_vec = [dV*cos(alpha)*v_minus_dir(1), dV*sin(alpha)*v_minus_dir(2), 0]'
            v_plus = dV_vec + obj.state(:,2)
            nstate = [obj.state(:,1), v_plus];
            
        
        
        end

        function obj = set.state(obj,new_state)

            obj.state = new_state;
        end

        



    end


end
