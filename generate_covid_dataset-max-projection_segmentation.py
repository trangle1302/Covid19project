import os
import hpacellseg.cellsegmentator as cellsegmentator
from hpacellseg.utils import label_cell, label_nuclei
import glob
import imageio
import numpy as np
from skimage import measure
from geojson import Polygon as geojson_polygon
from shapely.geometry import Polygon as shapely_polygon
from geojson import Feature, FeatureCollection, dump

def segment(image_folder, save_dir, segmentator):
    covid = glob.glob(image_folder + '/' + '*_red.tiff')
    er = [f.replace('red', 'yellow') for f in covid]
    nu = [f.replace('red', 'blue') for f in covid]
    images = [er, covid, nu]
    print(f'Found {len(images[2])} FOVs')
    img = imageio.imread(images[2][1])
    print(img.shape, img.max(), img.min())

    # For nuclei
    nuc_segmentations = segmentator.pred_nuclei(images[2])

    # For full cells
    cell_segmentations = segmentator.pred_cells(images)

    # post-processing
    for i, pred in enumerate(cell_segmentations):
        nuclei_mask, cell_mask = label_cell(nuc_segmentations[i], cell_segmentations[i])
        FOV_dir = os.path.join(save_dir, os.path.basename(covid[i]).replace("_red.tiff",""))
        if not os.path.exists(FOV_dir):
            os.makedirs(FOV_dir)
        print(nuclei_mask.dtype, cell_mask.max())
        #print(covid[i], os.path.join(FOV_dir,'dpnunet_cell_mask.png'))
        imageio.imwrite(os.path.join(FOV_dir,'dpnunet_cell_mask.png'), cell_mask) 

        #print(covid[i], os.path.join(FOV_dir,'dpnunet_nuclei_mask.png'))
        imageio.imwrite(os.path.join(FOV_dir,'dpnunet_nuclei_mask.png'), nuclei_mask) 

def masks_to_polygon(img_mask, label=None, save_name=None, simplify_tol=1.5):
    """
    example:
    label = 'cell'
    savename = os.path.join(sample_dest, "annotation.json")
    masks_to_polygon(img_mask, label=label, save_name=savename)
    """
    # for img_mask, for cells on border, should make sure on border pixels are # set to 0
    shape_x, shape_y = img_mask.shape
    shape_x, shape_y = shape_x - 1, shape_y - 1
    img_mask[0, :] = img_mask[:, 0] = img_mask[shape_x, :] = img_mask[
        :, shape_y
    ] = 0
    features = []
    # Get all object ids, remove 0 since this is background
    ind_objs = np.unique(img_mask)
    ind_objs = np.delete(ind_objs, np.where(ind_objs == 0))
    for obj_int in np.nditer(ind_objs):
        # Create binary mask for current object and find contour
        img_mask_loop = np.zeros((img_mask.shape[0], img_mask.shape[1]))
        img_mask_loop[img_mask == obj_int] = 1
        contours_find = measure.find_contours(img_mask_loop, 0.5)
        if len(contours_find) == 1:
            index = 0
        else:
            pixels = []
            for _, item in enumerate(contours_find):
                pixels.append(len(item))
            index = np.argmax(pixels)
        contour = contours_find[index]

        contour_as_numpy = contour[:, np.argsort([1, 0])]
        contour_as_numpy[:, 1] = np.array(
            [img_mask.shape[0] - h[0] for h in contour]
        )
        contour_asList = contour_as_numpy.tolist()

        if simplify_tol is not None:
            poly_shapely = shapely_polygon(contour_asList)
            poly_shapely_simple = poly_shapely.simplify(
                simplify_tol, preserve_topology=False
            )
            contour_asList = list(poly_shapely_simple.exterior.coords)
            contour_as_Numpy = np.asarray(contour_asList)

        # Create and append feature for geojson
        pol_loop = geojson_polygon([contour_asList])
        full_label = label + "_idx"
        index_number = int(obj_int - 1)
        features.append(
            Feature(
                geometry=pol_loop,
                properties={full_label: index_number, "label": label},
            )
        )

    if save_name:
        feature_collection = FeatureCollection(
            features, bbox=[0, 0, img_mask.shape[1] - 1, img_mask.shape[0] - 1]
        )
        with open(save_name, "w") as f:
            dump(feature_collection, f)
        os.chmod(save_name, 0o777)

if __name__ == "__main__":
    # generate csvs, including all-tasks.csv, all-well-tasks.csv
    # where to host the csvs
    #CSV_FOLDER = "/data/trang/covid19_data_CZ8746_max_projection/10"
    #SEGMENTATION_FOLDER = "/data/trang/covid19_data_CZ8746_annotation/segmentation2"
    #CHANNELS = ["_blue.tiff", "_red.tiff", "_green.tiff", "_yellow.tiff"]
    #IMG_FOLDER = "/data/trang/covid19_data_CZ8746_max_projection/10"
    ACQUIRED_DATA_PATH = "/data/trang/HPA_DV9903_Prescreen"
    IMG_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection/DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"
    SEGMENTATION_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection_annotation/segmentation/DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"
    NUC_MODEL = "./nuclei-model.pth"
    CELL_MODEL = "./cell-model.pth"
    segmentator = cellsegmentator.CellSegmentator(
        NUC_MODEL,
        CELL_MODEL,
        scale_factor=0.8,
        device="cuda",
        padding=False,
        multi_channel_model=True,
    )
    if not os.path.exists(SEGMENTATION_FOLDER):
        os.makedirs(SEGMENTATION_FOLDER)
    segment(IMG_FOLDER, SEGMENTATION_FOLDER, segmentator)
