 1488  sudo apt-get install ca-certificates curl gnupg
 2053  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/ --out_dir ../../output/mb --scan "manual" --mode "maxz" --make_preview --stitch_channels "ALEXA 488"
 2056  cd input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/
 2062  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_STRIATUM_POSITIVE_FULL_SLICE__2021-10-16T11_44_13-Measurement\ 1/ --out_dir ../../output/mb --scan "manual" --mode "maxz" --make_preview --stitch_channels "ALEXA 488"
 2066  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/Images/ --out_dir ../../output/mb --scan "manual" --mode "maxz" --make_preview --stitch_channels "ALEXA 488"
 2067  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/Images/ --out_dir ../../output/mb --scan "manual" --mode "maxz" --make_preview --stitch_channels "Alexa 488"
 2081  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/Images/ --out_dir ../../output/mb --scan "auto" --mode "maxz" --make_preview --stitch_channels "DAPI" "Alexa 488" "Alexa 647" --adaptive
 2084  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/Images/ --out_dir ../../output/mb --scan "auto" --mode "maxz" --make_preview --stitch_channels "DAPI" "Alexa 488" "Alexa 647" --adaptive
 2087  python -m biostitch --img_dir ../../input/egfp_images_230825/488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE__2021-10-16T10_14_08-Measurement\ 1/Images/ --out_dir ../../output/mb --scan "auto" --mode "maxz" --make_preview --stitch_channels "DAPI" "Alexa 488" "Alexa 647" --adaptive
 2092  ../../scripts/bftools/bfconvert -tilex 4096 -tiley 4096 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE.tif 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif
 2094  convert -resize 5% 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_resized_normal.tif
 2097  convert -resize 5% 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_resized_normal.tif
 2104  ../../scripts/bftools/showinf 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif 
 2106  ../../scripts/bftools/showinf -onexml-only 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif 
 2107  ../../scripts/bftools/showinf -omexml-only 488_TH_647_GFP_MIDBRAIN_POSITIVE_FULL_SLICE_normal.tif 
 2138  history | grep 488
 2143  history | grep 488 > phenix_imaging_commands.sh/
 2145  history | grep 488 > phenix_imaging_commands.sh
