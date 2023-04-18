import os
import imageio
import numpy as np
from glob import glob
import sys
import multiprocessing
import concurrent.futures
import pandas as pd
from matplotlib import pyplot as plt
from collections import OrderedDict
from skimage.filters import threshold_otsu

ACQUIRED_DATA_PATH = "/data/trang/HPA_DV9903_Prescreen"
PLATEID = "DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"
SAMPLES_DEST = f"{ACQUIRED_DATA_PATH}_max_projection_annotation/segmentation/{PLATEID}"
IMG_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection/{PLATEID}"
ALL_SAMPLES = []
for root, dirs, files in os.walk(SAMPLES_DEST):
    if ("dpnunet_nuclei_mask.png" in files) and (
        "dpnunet_cell_mask.png" in files
    ):
        img_path = os.path.join(IMG_FOLDER, os.path.basename(root))
        ALL_SAMPLES.append((img_path, root))

IMAGE_CHS = OrderedDict(
    {"protein": "_green.tiff", "virus": "_red.tiff", "er": "_yellow.tiff"}
)
MASK_FILES = OrderedDict(
    {
        "nuclei_mask": "dpnunet_nuclei_mask.png",
        "cell_mask": "dpnunet_cell_mask.png",
    }
)

CELL_IMAGE_INDEX = ["image_id", "cell_id"]
CHANNELS = ["protein", "virus", "er"]
AREAS_TO_QUANTIFY = ["nuclei", "cytosol", "cell"]
VALUES_TO_QUANTIFY = [
    "otsu_threshold",
    "integration",
    "integration_after_otsu_threshold",
    "mean",
    "mean_after_otsu_threshold",
    "median",
    "median_after_otsu_threshold",
]


def quantify_image(paths):
    image_path, mask_path = paths[0], paths[1]
    print("*" * 15, image_path, "*" * 15)
    try:
        image_channels = {}
        for key, value in IMAGE_CHS.items():
            image_channels[key] = imageio.imread(
                os.path.join(image_path + value)
            )

        masks = {}
        for key, value in MASK_FILES.items():
            masks[key] = imageio.imread(os.path.join(mask_path, value))
        masks["cytosol_mask"] = masks["cell_mask"] * (
            masks["cell_mask"] != masks["nuclei_mask"]
        )
        cell_indexes = np.unique(masks["nuclei_mask"])

        out = []
        if len(cell_indexes) < 2:
            with open(os.path.join(mask_path, "error3.txt"), "w") as f:
                f.write("no cell found")
            return

        for cell_idx in cell_indexes:
            if cell_idx == 0:
                continue
            current_cell_info = {}
            current_cell_info["image_id"] = image_path
            current_cell_info["cell_id"] = cell_idx

            for _, quantify_area in enumerate(
                AREAS_TO_QUANTIFY
            ):  # area ('nuclei', 'cytosol', 'cell')
                current_cell_area_mask_bool = (
                    masks[quantify_area + "_mask"] == cell_idx
                )
                current_cell_info[
                    "-".join([quantify_area, "size"])
                ] = np.count_nonzero(current_cell_area_mask_bool)
                # if the masked area does not exists, then abort this loop
                if not current_cell_info["-".join([quantify_area, "size"])]:
                    # out.append(current_cell_info)
                    with open(
                        os.path.join(mask_path, "error3.txt"), "w"
                    ) as f:
                        f.write("size zero, discarded")
                    return
                    # continue

                for _, channel in enumerate(
                    CHANNELS
                ):  # image channel ['protein', 'virus', 'er']
                    # only selected area remained
                    current_cell_area_masked = image_channels[channel][
                        current_cell_area_mask_bool
                    ]  # masked and flattend to one dim, only leave the area with values

                    # otsu_threshold_value = int(
                    #    threshold_otsu(current_cell_area_masked)
                    # )
                    # only values above otsu_threshold_value remained. one dim as well
                    # current_cell_area_otsued = current_cell_area_masked[
                    #    current_cell_area_masked > otsu_threshold_value
                    # ]

                    # current_cell_info[
                    #    "-".join([channel, quantify_area, "otsu_threshold"])
                    # ] = otsu_threshold_value
                    current_cell_info[
                        "-".join([channel, quantify_area, "integration"])
                    ] = np.sum(current_cell_area_masked)
                    # current_cell_info[
                    #    "-".join(
                    #        [
                    #            channel,
                    #            quantify_area,
                    #            "integration_after_otsu_threshold",
                    #        ]
                    #    )
                    # ] = np.sum(current_cell_area_otsued)
                    current_cell_info[
                        "-".join([channel, quantify_area, "min"])
                    ] = np.min(current_cell_area_masked).astype(int)
                    current_cell_info[
                        "-".join([channel, quantify_area, "max"])
                    ] = np.max(current_cell_area_masked).astype(int)
                    current_cell_info[
                        "-".join([channel, quantify_area, "std"])
                    ] = np.std(current_cell_area_masked).astype(int)
                    current_cell_info[
                        "-".join([channel, quantify_area, "mean"])
                    ] = np.mean(current_cell_area_masked).astype(int)
                    # get the mean value for the largest 600 pixel values
                    current_cell_info[
                        "-".join([channel, quantify_area, "largest500mean"])
                    ] = np.mean(
                        current_cell_area_masked[
                            np.argsort(current_cell_area_masked)[-600:]
                        ]
                    ).astype(
                        int
                    )
                    array_size = current_cell_area_masked.size
                    current_cell_info[
                        "-".join([channel, quantify_area, "50to75percentmean"])
                    ] = np.mean(
                        current_cell_area_masked[
                            np.argsort(current_cell_area_masked)[
                                int(array_size * 0.5) : int(array_size * 0.75)
                            ]
                        ]
                    ).astype(
                        int
                    )
                    # current_cell_info[
                    #    "-".join(
                    #        [
                    #            channel,
                    #            quantify_area,
                    #            "mean_after_otsu_threshold",
                    #        ]
                    #    )
                    # ] = np.mean(current_cell_area_otsued).astype(int)
                    current_cell_info[
                        "-".join([channel, quantify_area, "median"])
                    ] = np.median(current_cell_area_masked).astype(int)
                    # current_cell_info[
                    #    "-".join(
                    #        [
                    #            channel,
                    #            quantify_area,
                    #            "median_after_otsu_threshold",
                    #        ]
                    #    )
                    # ] = np.median(current_cell_area_otsued).astype(int)
            # append all values
            out.append(current_cell_info)
        df = pd.DataFrame(out)
        if os.path.exists(os.path.join(mask_path, "cell_quantify.csv")):
            os.remove(os.path.join(mask_path, "cell_quantify.csv"))
        csv_file = os.path.join(mask_path, "cell_quantify.csv")
        df.to_csv(csv_file, index=False)
        # SAMPLES_CSV.loc[
        #    SAMPLES_CSV["sample"] == image_path, "cell_quantified"
        # ] = True
        # SAMPLES_CSV.to_csv(
        #    os.path.join(SAMPLES_DEST, "tasks.csv"), index=False
        # )
    except Exception as e:
        print("$" * 20, "error: ", mask_path, str(e))
        with open(os.path.join(mask_path, "error.txt"), "w") as f:
            f.write(str(e))
        # pass
        # SAMPLES_CSV.loc[
        #    SAMPLES_CSV["sample"] == image_path, "cell_quantified"
        # ] = "Error"
        # SAMPLES_CSV.to_csv(
        #    os.path.join(SAMPLES_DEST, "tasks.csv"), index=False
        # )


def quantify_all():
    # same the following for next round
    # for sample in ALL_SAMPLES:
    #    if os.path.exists(os.path.join(sample, "cell_quantify.csv")):
    #        SAMPLES_CSV.loc[
    #            SAMPLES_CSV["sample"] == sample, "cell_quantified"
    #        ] = 1
    #    if os.path.exists(os.path.join(sample, "error2.txt")):
    #        SAMPLES_CSV.loc[
    #            SAMPLES_CSV["sample"] == sample, "cell_quantified"
    #        ] = 0
    # SAMPLES_CSV.to_csv(os.path.join(SAMPLES_DEST, "tasks.csv"), index=False)
    with multiprocessing.Pool() as pool:
        pool.map(quantify_image, ALL_SAMPLES)
    # for sample in ALL_SAMPLES:
    #    quantify_image(sample)


if __name__ == "__main__":
    print(ALL_SAMPLES[:9])
    quantify_all()
