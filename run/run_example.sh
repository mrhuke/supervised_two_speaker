readonly command="addpath('./snr'); addpath('./cochleagram'); sig=load('sample/t11_lwby6p_m30_lrwp2a.-9dB.val2'); rmask=twoSpeaker(sig, [11,30], 'acoustDym_iter', 256, 16, .5, 3, '.');"
matlab -nodesktop -nosplash -nodisplay -r "$command"
