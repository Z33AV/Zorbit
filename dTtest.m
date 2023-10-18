Re = 6378.14;
orb = conic(4.5*Re,0.75,0,0,0,'earth');
sc = spacecraft_z(orb,[],'XYZ','earth');
sc.TA = deg2rad(90);
sc.E = spacecraft_z.TA2E(sc.TA,sc.orbit.e)
TA_apo = deg2rad(180);

sc.time_till(TA_apo)