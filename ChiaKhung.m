function FrameMat = ChiaKhung(WavFile,FrameSize)
%H�M N�Y D�NG ?? CHIA WAVFILE BAN ??U TH�NH C�C KHUNG
%--------------------------------------------------------------------------------------------------
%Tham s? ??u v�o: +) WavFile: gi� tr? bi�n ?? c?a c�c ph?n t? trong wavfile ban ??u
%                 +) FrameSize: S? l??ng ph?n t? c?a 1 khung
%Tham s? ??u ra: FrameMat: M?ng ch?a c�c khung trong ?�
%		                 + M?i c?t l� 1 khung
%                        + M?i h�ng trong 1 c?t t??ng ?ng l� c�c ph?n t? c?a khung
%--------------------------------------------------------------------------------------------------
%B??C 1: ?? D�I WAVFILE L�M �TR�N� L�N S? L?N H?N (HO?C B?NG), G?N NH?T, CHIA H?T CHO FRAMESIZE
SoDu = mod(length(WavFile),FrameSize);
if(SoDu > 0)
        NewlengthwavFile = length(WavFile)+FrameSize-SoDu;  %NewlengthwavFile ?? l?u s? v?a l�m �tr�n�
else 
        NewlengthwavFile = length(WavFile);
end
%--------------------------------------------------------------------------------------------------
% B??C 2: D�NG 1 M?NG TRUNG GIAN (C� NewlengthwavFile PH?N T?) CH?A M?NG WAVFILE
%(CH?A GI� TR? C�C BI�N ??) V� CH�N TH�M C�C S? 0 V�O TR??C HO?C SAU 
SoLuongKhung = 2*NewlengthwavFile/FrameSize;
%N?u s? l??ng ph?n t? c?n th�m >= ( frameSize/2)
if(FrameSize-SoDu >= FrameSize/2 )
    % M?ng trung gian d�ng ?? l?p ch?a WavFile v� l?p ??y b?ng nh?ng gi� tr? 0 
    TrungGian = [zeros(FrameSize/2,1);WavFile;zeros(FrameSize/2-SoDu,1)];   
else     %N?u s? l??ng ph?n t? c?n th�m < ( frameSize/2)
    TrungGian = [zeros(FrameSize/2,1);WavFile;zeros(FrameSize-SoDu,1)];
end
% T?o M?ng fameMat
k = 0; 
for i = 1:2:SoLuongKhung
    for j = 1:FrameSize
        % Gi� tr? frameSize*k ?? kh?i t?o gi� tr? b?t ??u cho 1 khung m?i c?a khung l?
        FrameMat(j,i) = TrungGian(FrameSize*k+j);    
        if(FrameSize*k+j+FrameSize/2 <= NewlengthwavFile)
            % Gi� tr? frameSize*k+frameSize/2 ?? kh?i t?o gi� tr? b?t ??u 1 khung m?i c?a khung ch?n
            FrameMat(j,i+1) = TrungGian(FrameSize*k+FrameSize/2+j);    
        end
    end
    k = k+1;
end

