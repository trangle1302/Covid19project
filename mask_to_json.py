import os
import numpy as np
from skimage import measure
from PIL import Image
from geojson import Polygon as geojson_polygon
from shapely.geometry import Polygon as shapely_polygon
from geojson import Feature, FeatureCollection, dump
import pandas as pd
from shutil import copyfile, rmtree
from imageio import imread


SEGMENTATION_FOLDER = "/data/trang/covid19_data_CZ8746_annotation/segmentation2"

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
    
    well_folders = os.listdir(SEGMENTATION_FOLDER)
    label = 'cell'
    i = 0
    for well_folder in well_folders:
        if os.path.exists(os.path.join(SEGMENTATION_FOLDER, well_folder, 'dpnunet_cell_mask.png')):
            try:
                #print('*'*10, root)
                img_mask = imread(os.path.join(SEGMENTATION_FOLDER, well_folder, 'dpnunet_cell_mask.png'))
                save_name = os.path.join(os.path.join(SEGMENTATION_FOLDER, well_folder, "annotation.json"))
                
                masks_to_polygon(img_mask, label, save_name)
            except:
                i += 1
                print(save_name)
    print(i)
