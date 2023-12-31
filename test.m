clear all
clc
Re = 6378.14;


%Validating against AAE 532 PS6 PRoblem 3
r = [7*Re; 2*Re; 3*Re];
v = [3.4; -0.2; 0.1];

state = [r,v];

sc = spacecraft_z(0, state, 'xyz','earth')

[a,e,w,O,i,E,h,gamma,M] = sc.kepels()
sc.E = E;
sc.h = h;
sc.gamma = gamma
sc.M = M;
sc.orbit = conic(a,e,w,O,i,sc.body)
sc.TA = spacecraft_z.E2TA(sc.E,sc.orbit.e)

TA_targ = deg2rad(234.5);
dT = sc.time_till(TA_targ)

sc.orbit.plot2(0,360,1);
sc.orbit.plot3(0,360,0,2);

nstate = sc.impulse(0.3,deg2rad(65));

sc.state = nstate;
[a,e,w,O,i,E,h,gamma,M] = sc.kepels()
sc.E = E;
sc.h = h;
sc.gamma = gamma
sc.M = M;
sc.orbit = conic(a,e,w,O,i,sc.body)
sc.TA = spacecraft_z.E2TA(sc.E,sc.orbit.e)
sc.orbit.plot2(0,360,1);
sc.orbit.plot3(0,360,0,2);
