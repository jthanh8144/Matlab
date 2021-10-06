% Nguong la Hz
function [F0, x, yy] = AMDF(y, Fs, khung, hzMin, hzMax, nguong)
    timeMin = floor(Fs/hzMax); %
    timeMax = floor(Fs/hzMin); 
    khungMau = khung*Fs;
    index = 1;
    F0 = zeros(1,floor(length(y)/khungMau));
    F0 = 1/Fs;
    i = 0;
    while index < length(y)
        if(index + khungMau < length(y))
            y1 = y(index:index + khungMau);

            yy = zeros(1,2*length(y1)); %Hàm tu tuong quan
            for k = -length(y1):length(y1)
                for j = 1:length(y1)
                    if (j - k > 1) && (j - k < length(y1) - 1)
                        %Cong don tich tin hieu y[j]*y[j - k]
                        %Voi k la so mau bi lech
                        yy(k + length(y1) + 1) = yy(k + length(y1) + 1) + abs(y1(j) - y1(j - k));
                    end
                end
            end
            i = i + 1;
            %Tìm F0 max
            valueMin = max(yy);
            time = 0;
            for j = length(yy)/2 + timeMin:length(yy)/2 + timeMax
                if (yy(j) < valueMin) && yy(j) < max(yy)*nguong
                    valueMin = yy(j);
                    time = j - length(yy)/2;
                end
            end
            F0(i) = Fs/time;
        end
        index = index + khungMau/2 + 1;
    end
    x = [1:length(F0)];
    x = (x./Fs)*(length(y)/length(F0));
end