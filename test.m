clear all
clc
Re = 6378.14;


%Validating against AAE 532 PS6 PRoblem 3
r = [7*Re; 2*Re; 3*Re];
v = [3.4; -0.2; 0.1];

state = [r,v];

sc = spacecraft_z(0, state, 'xyz','earth')

[a,e,w,O,i,E,h,gamma,M] = sc.kepels()
sc.orbit = conic(a,e,w,O,i,sc.body)

sc.orbit.plot2(0,360);
sc.orbit.plot3(0,360,0);
