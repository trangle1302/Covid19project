from collections import OrderedDict


EXPERIMENT = "CZ8780"
ACQUIRED_DATA_PATH = "/data/trang/HPA_CZ8780"
PLATEID = "CZ8780_plate_II"
FORMATTED_DATA_PATH = f"{ACQUIRED_DATA_PATH}_formatted"
ANNOTATION_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection_annotation"
SAMPLES_DEST = f"{ACQUIRED_DATA_PATH}_max_projection_annotation/segmentation/{PLATEID}"
IMG_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection/{PLATEID}"
SEGMENTATION_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection_annotation/segmentation/CZ8780_plate_I"#DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"

NUC_MODEL = "./nuclei-model.pth"
CELL_MODEL = "./cell-model.pth"


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
THRESHOLD = False # whether of otsu threshold to reduce noise