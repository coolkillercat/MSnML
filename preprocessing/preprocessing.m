datasets = {'G2S2(3)', 'G2FS2(3)','G2S1(6)', 'G2S2(6)', 'G2FS1(6)','G2FS2(6)','G2FS2(6)_REDO','ControlFetuin','ControlTransferrin','Alpha.2,3Neur.Fetuin'};
ratio = [1,1,0,0,0,0,0,-1,-1, -1];
dmass = CMass.get_atom_mass('Li_fix');
MS_level = 'MS5';
path = 'D:\TPP\data';
%
if strcmp(MS_level,'MS5') == 1
    resolution = 10000;
else
    resolution = 23000;
end
for di = 1:length(datasets)
    dataset = datasets{di};
    path = [path,'\',dataset,'\',MS_level,'\'];
    MSfiles = dir([path, '*.mzXML']);
    MSns = {};
    if ratio(di) == 1
        load('template_theoretical_align_3.mat')
    elseif ratio(di) == 0
        load('template_theoretical_align_6.mat')
    else
        load('template_mix.mat')
    end
    for i = 1:length(MSfiles)
        filename = [path, MSfiles(i).name];
        MSns{end+1} = CMSn(filename,'resolution',resolution,'template',template, 'ratio', ratio(di));
    end
    save([dataset,'_',MS_level,'.mat'],'MSns')
end
glycan_list = datasets;
MS_level = 'MS5';
for gi = 1:length(glycan_list)
    IntensityMatrix = [];
    load([glycan_list{gi}, '_',MS_level,'.mat'])
    ratio3 = [];
    for i = 1:length(MSns)
        ratio3(end+1) = MSns{i}.mRatio;
        if max(MSns{i}.mSumIntensity) == 0 
            IntensityMatrix = [IntensityMatrix zeros(13333,1)];
            continue;
        end
        IntensityMatrix = [IntensityMatrix MSns{i}.mSumIntensity/max(MSns{i}.mSumIntensity)];
    end
    mass_features = MSns{1}.mM;
    
window = string.empty;
method = string.empty;
energy = [];
rr = string.empty;
for i = 1:length(MSns)
    ms5 = MSns{i};
    window(end+1) = ms5.mType;
    method(end+1) = ms5.mExperiment;
    energy(end+1) = ms5.mCollisionEnergy;
    rr(end+1) = convertCharsToStrings(ms5.mMixRatio);
end
metadata = table(window', method',energy', rr');

    folder_name = ['test_', glycan_list{gi}, '_',MS_level,'\'];
    mkdir(folder_name);
    savepath = [folder_name];
    writetable(metadata, [savepath,'metadata.csv'])
    writematrix(IntensityMatrix, [savepath, 'input.csv']);
    writematrix(ratio3, [savepath,'output.csv']);
    writematrix(mass_features, [savepath,'mass_features.csv'])
end