%% K = 4
i = 1;
filterbankSettings{i}.numberOfChannels = 4;
filterbankSettings{i}.decimationFactor = 2;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 16;
filterbankSettings{i}.latency = 579;
filterbankSettings{i}.peaqScore = -0.35;
filterbankSettings{i}.signalToAliasRatio = 55;

i = i+1;
filterbankSettings{i}.numberOfChannels = 4;
filterbankSettings{i}.decimationFactor = 3;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 32;
filterbankSettings{i}.latency = 593;
filterbankSettings{i}.peaqScore = -1.53;
filterbankSettings{i}.signalToAliasRatio = 38.5;

i = i+1;
filterbankSettings{i}.numberOfChannels = 4;
filterbankSettings{i}.decimationFactor = 3;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 64;
filterbankSettings{i}.latency = 625;
filterbankSettings{i}.peaqScore = 0.07;
filterbankSettings{i}.signalToAliasRatio = 62;

%% K = 8
i = i+1;
filterbankSettings{i}.numberOfChannels = 8;
filterbankSettings{i}.decimationFactor = 4;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 32;
filterbankSettings{i}.latency = 593;
filterbankSettings{i}.peaqScore = -2.33;
filterbankSettings{i}.signalToAliasRatio = 53;

i = i+1;
filterbankSettings{i}.numberOfChannels = 8;
filterbankSettings{i}.decimationFactor = 4;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 64;
filterbankSettings{i}.latency = 625;
filterbankSettings{i}.peaqScore = 0.05;
filterbankSettings{i}.signalToAliasRatio = 95;

i = i+1;
filterbankSettings{i}.numberOfChannels = 8;
filterbankSettings{i}.decimationFactor = 5;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 64;
filterbankSettings{i}.latency = 621;
filterbankSettings{i}.peaqScore = 3.58e-4;
filterbankSettings{i}.signalToAliasRatio = 64;

i = i+1;
filterbankSettings{i}.numberOfChannels = 8;
filterbankSettings{i}.decimationFactor = 6;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 128;
filterbankSettings{i}.latency = 683;
filterbankSettings{i}.peaqScore = 0.02;
filterbankSettings{i}.signalToAliasRatio = 78;

%% K = 16
i = i+1;
filterbankSettings{i}.numberOfChannels = 16;
filterbankSettings{i}.decimationFactor = 8;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 64;
filterbankSettings{i}.latency = 617;
filterbankSettings{i}.peaqScore = -1.65;
filterbankSettings{i}.signalToAliasRatio = 52;

i = i+1;
filterbankSettings{i}.numberOfChannels = 16;
filterbankSettings{i}.decimationFactor = 11;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 128;
filterbankSettings{i}.latency = 677;
filterbankSettings{i}.peaqScore = -0.58;
filterbankSettings{i}.signalToAliasRatio = 50;

i = i+1;
filterbankSettings{i}.numberOfChannels = 16;
filterbankSettings{i}.decimationFactor = 13;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 256;
filterbankSettings{i}.latency = 803;
filterbankSettings{i}.peaqScore = -0.72;
filterbankSettings{i}.signalToAliasRatio = 55;

%% K = 32
i = i+1;
filterbankSettings{i}.numberOfChannels = 32;
filterbankSettings{i}.decimationFactor = 16;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 128;
filterbankSettings{i}.latency = 673;
filterbankSettings{i}.peaqScore = -1.47;
filterbankSettings{i}.signalToAliasRatio = 49;

i = i+1;
filterbankSettings{i}.numberOfChannels = 32;
filterbankSettings{i}.decimationFactor = 21;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 256;
filterbankSettings{i}.latency = 799;
filterbankSettings{i}.peaqScore = -0.26;
filterbankSettings{i}.signalToAliasRatio = 53;

i = i+1;
filterbankSettings{i}.numberOfChannels = 32;
filterbankSettings{i}.decimationFactor = 26;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 512;
filterbankSettings{i}.latency = 1033;
filterbankSettings{i}.peaqScore = -0.17;
filterbankSettings{i}.signalToAliasRatio = 58;

%% K = 64
i = i+1;
filterbankSettings{i}.numberOfChannels = 64;
filterbankSettings{i}.decimationFactor = 32;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 256;
filterbankSettings{i}.latency = 769;
filterbankSettings{i}.peaqScore = -1.10;
filterbankSettings{i}.signalToAliasRatio = 48;

i = i+1;
filterbankSettings{i}.numberOfChannels = 64;
filterbankSettings{i}.decimationFactor = 43;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 512;
filterbankSettings{i}.latency = 1015;
filterbankSettings{i}.peaqScore = -0.34;
filterbankSettings{i}.signalToAliasRatio = 52;

i = i+1;
filterbankSettings{i}.numberOfChannels = 64;
filterbankSettings{i}.decimationFactor = 51;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 1024;
filterbankSettings{i}.latency = 1492;
filterbankSettings{i}.peaqScore = -0.47;
filterbankSettings{i}.signalToAliasRatio = 55;

%% K = 128
i = i+1;
filterbankSettings{i}.numberOfChannels = 128;
filterbankSettings{i}.decimationFactor = 64;
filterbankSettings{i}.oversamplingRate = 2;
filterbankSettings{i}.filterLength = 512;
filterbankSettings{i}.latency = 993;
filterbankSettings{i}.peaqScore = -0.75;
filterbankSettings{i}.signalToAliasRatio = 47;

i = i+1;
filterbankSettings{i}.numberOfChannels = 128;
filterbankSettings{i}.decimationFactor = 85;
filterbankSettings{i}.oversamplingRate = 3/2;
filterbankSettings{i}.filterLength = 1024;
filterbankSettings{i}.latency = 1441;
filterbankSettings{i}.peaqScore = -0.35;
filterbankSettings{i}.signalToAliasRatio = 50;

i = i+1;
filterbankSettings{i}.numberOfChannels = 128;
filterbankSettings{i}.decimationFactor = 102;
filterbankSettings{i}.oversamplingRate = 5/4;
filterbankSettings{i}.filterLength = 2048;
filterbankSettings{i}.latency = 2465;
filterbankSettings{i}.peaqScore = -0.41;
filterbankSettings{i}.signalToAliasRatio = 63;