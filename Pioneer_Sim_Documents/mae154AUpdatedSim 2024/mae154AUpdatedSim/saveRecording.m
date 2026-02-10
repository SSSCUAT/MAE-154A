load dataRecorder.mat

now = clock;
fileName = ['simRecording_' num2str(now(1)) num2str(now(2),'%02.f') num2str(now(3),'%02.f') ...
    '_' num2str(now(4),'%02.f') num2str(now(5),'%02.f') '.txt'];

fileID = fopen(fileName,'w');

for i=1:length(recordedData)
    for j=1:16
        if j<16
            fprintf(fileID,'%5d ',recordedData(j,i));
        else
            fprintf(fileID,'%5d \n',recordedData(j,i));
        end
%            fprintf(fileID,'%5d \n',recordedData(j,i));
    end
    
end
fclose(fileID);
