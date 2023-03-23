#!/bin/bash -l
ffmpeg -framerate 4 -i gganim_plot%04d.png -vcodec libx264 -pix_fmt yuv420p pheno_mean.mp4
