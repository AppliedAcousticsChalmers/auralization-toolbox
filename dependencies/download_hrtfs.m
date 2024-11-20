function [] = download_hrtfs(hrir_path)
% download HRIRs (skipped automatically if the HRIR dataset already exists)

if ~isfile(hrir_path)
    hrir_path_components = strsplit(hrir_path, '/');
    for n = 1 : length(hrir_path_components)-1
        mkdir(hrir_path_components{n});
    end
    fprintf('Cannot find the HRTFs. '); 
    fprintf('Downloading them from https://zenodo.org/record/3928297 ... ');
    websave(hrir_path, 'https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1');
    fprintf('done.\n\n');
end

end