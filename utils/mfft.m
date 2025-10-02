function f=mfft(data,dir)


    if dir ==1
        f=ifftshift(ifft(fftshift(data,1),[],1),1);
    elseif dir==2
        f=ifftshift(ifft(fftshift(data,2),[],2),2);
    else
        f=ifftshift(ifft2(fftshift(data)));
    end

end
