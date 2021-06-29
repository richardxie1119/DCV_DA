function [fids,calib] = loadBrukerFIDs(data_dir,read_size,fid_index,fid_length)
    fids = [];
    calib = [];
    file_dir = dir(strcat(data_dir,'/**/*.d'));
    file_num = size(file_dir);
    for i = 1:file_num(1)
        f_id = fopen(strcat(file_dir(i).folder,'/',file_dir(i).name,'/ser'));
        for j = 1:length(fid_index)
            disp(strcat('loading FID...',int2str(j),'-',int2str(fid_index(j))))
            fseek(f_id,4*(fid_index(j)-1)*fid_length,'bof');
            fid = fread(f_id,fid_length,'int32');

            if read_size == 'all'
                fids = [fids,fid];
            else
                fids = [fids,fid(1:read_size)];
            end
        end
        
        fclose(f_id);
       
        method_dir = dir(strcat(file_dir(i).folder,'/',file_dir(i).name,'/*.m'));
        method = fileread(strcat(method_dir(1).folder,'/',method_dir(1).name,'/apexAcquisition.method'));
        matchStr = regexp(method,'<param name="ML1"><value>(.+?)</value>','tokens');
        calibA = matchStr{1};
        matchStr2 = regexp(method,'<param name="ML2"><value>(.+?)</value>','tokens');
        calibB = matchStr2{1};
        matchStr3 = regexp(method,'<param name="ML3"><value>(.+?)</value>','tokens');
        calibC = matchStr3{1};

        calib = [calib;str2double(calibA),str2double(calibB),str2double(calibC)];
    end
    
    fids = fids.';
    
end

