
prm.defaultFrameLength = 4.0;
prm.F0frameUpdateInterval=5.000000;
prm.F0archUpperBound=330;
prm.F0archLowerBound=50;
prm.spectralUpdateInterval=5.000000;

[x,fs]=audioread('jal_in_01_3.wav');
[f0, ap] = exstraightsource(x,fs, prm);
[sp] = exstraightspec(x, f0, fs, prm);


sy = exstraightsynth(f0, sp, ap, fs, prm);

audiowrite( 'test_dct_39.wav', sy/max(abs(sy)), 16000 );
