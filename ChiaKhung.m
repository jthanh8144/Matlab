function FrameMat = ChiaKhung(WavFile,FrameSize)
%HÀM NÀY DÙNG ?? CHIA WAVFILE BAN ??U THÀNH CÁC KHUNG
%--------------------------------------------------------------------------------------------------
%Tham s? ??u vào: +) WavFile: giá tr? biên ?? c?a các ph?n t? trong wavfile ban ??u
%                 +) FrameSize: S? l??ng ph?n t? c?a 1 khung
%Tham s? ??u ra: FrameMat: M?ng ch?a các khung trong ?ó
%		                 + M?i c?t là 1 khung
%                        + M?i hàng trong 1 c?t t??ng ?ng là các ph?n t? c?a khung
%--------------------------------------------------------------------------------------------------
%B??C 1: ?? DÀI WAVFILE LÀM “TRÒN” LÊN S? L?N H?N (HO?C B?NG), G?N NH?T, CHIA H?T CHO FRAMESIZE
SoDu = mod(length(WavFile),FrameSize);
if(SoDu > 0)
        NewlengthwavFile = length(WavFile)+FrameSize-SoDu;  %NewlengthwavFile ?? l?u s? v?a làm “tròn”
else 
        NewlengthwavFile = length(WavFile);
end
%--------------------------------------------------------------------------------------------------
% B??C 2: DÙNG 1 M?NG TRUNG GIAN (CÓ NewlengthwavFile PH?N T?) CH?A M?NG WAVFILE
%(CH?A GIÁ TR? CÁC BIÊN ??) VÀ CHÈN THÊM CÁC S? 0 VÀO TR??C HO?C SAU 
SoLuongKhung = 2*NewlengthwavFile/FrameSize;
%N?u s? l??ng ph?n t? c?n thêm >= ( frameSize/2)
if(FrameSize-SoDu >= FrameSize/2 )
    % M?ng trung gian dùng ?? l?p ch?a WavFile và l?p ??y b?ng nh?ng giá tr? 0 
    TrungGian = [zeros(FrameSize/2,1);WavFile;zeros(FrameSize/2-SoDu,1)];   
else     %N?u s? l??ng ph?n t? c?n thêm < ( frameSize/2)
    TrungGian = [zeros(FrameSize/2,1);WavFile;zeros(FrameSize-SoDu,1)];
end
% T?o M?ng fameMat
k = 0; 
for i = 1:2:SoLuongKhung
    for j = 1:FrameSize
        % Giá tr? frameSize*k ?? kh?i t?o giá tr? b?t ??u cho 1 khung m?i c?a khung l?
        FrameMat(j,i) = TrungGian(FrameSize*k+j);    
        if(FrameSize*k+j+FrameSize/2 <= NewlengthwavFile)
            % Giá tr? frameSize*k+frameSize/2 ?? kh?i t?o giá tr? b?t ??u 1 khung m?i c?a khung ch?n
            FrameMat(j,i+1) = TrungGian(FrameSize*k+FrameSize/2+j);    
        end
    end
    k = k+1;
end

