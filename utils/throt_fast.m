function Rth = throt_fast(phi, theta)
    % Rotation um Achse im Azimut theta
    Rth = zrot(theta) * xrot(phi) * zrot(-theta);
end
