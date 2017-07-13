function [idx0, idx1, nr1, riv_msk] = riv_mask(riv_msk_fil)
riv_msk = edge(riv_msk);
idx1 = find(riv_msk_fil == 1); idx0 = find(riv_msk_fil == 0); nr1 = length(idx1);
