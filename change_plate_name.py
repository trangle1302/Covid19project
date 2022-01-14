import os
import glob

old_name = "1"
new_name = "10"

"""
SAMPLES_DIR = "/data/trang/covid19_data_CZ8746_max_projection/10/"
imlist = glob.glob(SAMPLES_DIR + "*.tiff")
print(len(imlist))
for f in imlist:
    name0 = os.path.basename(f)
    plate, well, sample, ch = name0.split("_")
    assert(plate == old_name)
    name1 = os.path.join(SAMPLES_DIR, "_".join((new_name, well, sample, ch)))
    print(f, name1)
    os.rename(f, name1)


SAMPLES_DEST = "/data/trang/covid19_data_CZ8746_annotation/segmentation"
paths = os.listdir(SAMPLES_DEST)
for d in paths:
    name0 = os.path.join(SAMPLES_DEST,d)
    plate, well, sample = d.split("_")
    assert(plate == old_name)
    name1 = os.path.join(SAMPLES_DEST, "_".join((new_name, well, sample)))
    print(name0, name1)
    os.rename(name0, name1)


# switch green and yellow
SAMPLES_DIR = "/data/trang/covid19_data_CZ8746_max_projection/10/"
imlist = glob.glob(SAMPLES_DIR + "*_yellow.tiff")
print(len(imlist))
for ori_yellow in imlist:
    tmp = ori_yellow.replace("_yellow.tiff", "_tmp.tiff")
    print(f"Moving {ori_yellow}.......{tmp}")
    os.rename(ori_yellow, tmp)

    ori_green = ori_yellow.replace("_yellow.tiff", "_green.tiff")
    print(f"Moving {ori_green}.......{ori_yellow}")
    os.rename(ori_green, ori_yellow)
    
    print(f"Moving {tmp}.......{ori_green}")
    os.rename(tmp, ori_green)

"""