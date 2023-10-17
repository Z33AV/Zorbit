clear all
close all
clc

%givens
Re = 6378.14
a = 1.4446*Re;
e = .22571;
i = 0.3;
O = 0.5;
w = 0.5;

init_orb = conic(a,e,w,O,i,'earth');
rmag_init = 1.65*Re;
vmag_init = 5.7;
FPA_init = deg2rad(-10.2); %rad

%generating state vector
r_init_rth =  [rmag_init,0,0];
v_init_rth = [vmag_init*sin(FPA_init), vmag_init*cos(FPA_init), 0]
p = a*(1-e^2);
theta = acos( p/rmag_init - 1 ) + w;
xyz_rth = angle2dcm(O,i,theta, 'ZXZ');
rth_xyz = xyz_rth';

r_init_xyz = rth_xyz*r_init_rth';
v_init_xyz = rth_xyz*v_init_rth';

%initializing spacecraft object

state = [r_init_xyz, v_init_xyz]
sc = spacecraft_z(init_orb,state,'xyz','earth')


%gut check plot initial orbit
sc.orbit.plot3(0,360,0);

%impulse test
nstate = sc.impulse(1.2,deg2rad(60));
sc.state = nstate
[a,e,w,O,i,E,h,gamma,M] = sc.kepels();
sc.orbit  = conic(a,e,w,O,i,'earth');
sc.orbit.plot2(0,360);
sc.orbit.plot3(0,360,0);









