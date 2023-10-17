clear all
close all
clc

orbit = conic(24555.824,0.61038,0,0,0,'earth');

sctest = spacecraft_z(orbit, [-29998.073,15394.196; -2.322,-1.421], 'eph','earth')
sctest.E = deg2rad(127.679);
dt = 5*60*60;
E1 = deg2rad(239.179);

fgmat = sctest.fgcalc(E1,dt);