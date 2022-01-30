# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# ### This is to format the acquired image plates to the format of covid-annotation and then generate the all-tasks.csv and all-well-tasks.csv
# ##### the csv files should be copied to the covid-annotation folder
import os
from shutil import copyfile
from glob import glob
import numpy as np
from imageio import imread, imsave


def format_covid_data(acquired_data_path, formatted_data_folder):
    """not fov_index_length should be set depending on how many fovs in one well. like 9 fov and 3 z-slice image, 3*9=27, then fov_index_length = 2. this is to ensure the z-stack be together in the annotation tool"""
    os.makedirs(formatted_data_folder, exist_ok=True)
    plates = list(
        filter(
            lambda item: item if not item.startswith(".") else None,
            os.listdir(acquired_data_path),
        )
    )
    row_mapping = {
        "r01": "A",
        "r02": "B",
        "r03": "C",
        "r04": "D",
        "r05": "E",
        "r06": "F",
        "r07": "G",
        "r08": "H",
    }
    channel_mapping = {
        "ch1": "blue.tiff", # Hoechst
        "ch2": "green.tiff", # Alexa 488 - protein
        "ch3": "red.tiff", # Alexa 555 - covid
        "ch4": "yellow.tiff", # Alexa 647 - ER
    }
    for plate in plates:
        plate_folder = os.path.join(formatted_data_folder, plate)
        os.makedirs(plate_folder, exist_ok=True)
        all_images = []
        plate_path = os.path.join(acquired_data_path, plate)
        image_folder = os.path.join(plate_path, "Images")
        print(image_folder)
        assert os.path.exists(image_folder)
        all_images = list(
            filter(
                lambda item: item if item.endswith(".tiff") else None,
                os.listdir(image_folder),
            )
        )
        fovs = list(map(lambda item: item.split("-")[0], all_images))
        fovs = list(set(fovs))
        fovs.sort()
        mapping_dict = {}
        for i, item in enumerate(fovs):
            row_index = row_mapping[item[:3]]
            column_index = str(int(item[4:6]))
            fov_index = str(int(item[7:9]))
            plane_index = str(int(item[10:12]))
            mapping_dict[item] = "".join(
                [row_index, column_index, "_", fov_index, "_", plane_index]
            )
            print(''.join([row_index, column_index, '_', str(fov_index)]))
        for image in all_images:
            source_path = os.path.join(
                acquired_data_path, plate, "Images", image
            )
            image_items = image.split("-")
            assert len(image_items) == 2 and (
                image_items[1][:3] in ["ch1", "ch2", "ch3", "ch4"]
            )
            fov_index = mapping_dict[image_items[0]]
            channel_index = channel_mapping[image_items[1][:3]]
            file_name = "_".join([plate, fov_index, channel_index])
            # print("{}:{}".format(image, file_name))
            dest_path = os.path.join(plate_folder, file_name)
            print("{}::::{}".format(source_path, dest_path))
            copyfile(source_path, dest_path)
    # print(len(all_images))


def generate_task_lists(csv_folder, data_base):
    images_found = []
    for root, dirs, files in os.walk(data_base):
        for file in files:
            if file.startswith(".") or not file.endswith(".tiff"):
                continue
            else:
                images_found.append(file)
    all_tasks = []
    all_well_tasks = []
    for image in images_found:
        image_items = image.split("_")
        all_tasks.append("_".join(image_items[:3]))
        all_well_tasks.append("_".join(image_items[:2]))
    all_tasks = list(set(all_tasks))
    all_well_tasks = list(set(all_well_tasks))
    all_tasks.sort()
    all_well_tasks.sort()
    well_tasks_file = open(os.path.join(csv_folder, "all-well-tasks.csv"), "w")
    all_tasks_file = open(os.path.join(csv_folder, "all-tasks.csv"), "w")
    well_tasks_file.write("task_id" + "\n")
    all_tasks_file.write("task_id" + "\n")
    for task in all_well_tasks:
        well_tasks_file.write(task + "\n")
    well_tasks_file.close()
    for task in all_tasks:
        all_tasks_file.write(task + "\n")
    all_tasks_file.close()


def proj_channel(channel_imgs, img_dtype):
    channel_imgs = np.max(
        [imread(channel_img) for channel_img in channel_imgs], axis=0
    )
    channel_imgs = np.asarray(channel_imgs, dtype=img_dtype)
    return channel_imgs


def generate_max_projection_images(formatted_data_folder, data_base):
    os.makedirs(data_base, exist_ok=True)
    plates = list(
        filter(
            lambda item: item if not item.startswith(".") else None,
            os.listdir(formatted_data_folder),
        )
    )
    for plate in plates:
        plate_folder = os.path.join(data_base, plate)
        os.makedirs(plate_folder, exist_ok=True)
        all_images = []
        image_folder = os.path.join(formatted_data_folder, plate)
        assert os.path.exists(image_folder)
        all_images = list(
            filter(
                lambda item: item if item.endswith(channels[0]) else None,
                os.listdir(image_folder),
            )
        )
        img_sample = imread(os.path.join(image_folder, all_images[0]))
        img_dtype = img_sample.dtype
        print(img_sample.shape)
        # print('&'*80)
        # print(img_dtype)
        fovs = list(
            map(lambda item: "_".join(item.split("_")[:3]), all_images)
        )
        fovs = list(set(fovs))
        # print(fovs)
        for _, fov in enumerate(fovs):
            fov_chs = []
            for i in range(len(channels)):
                fov_chs.append([])
                fov_chs[i] = glob(image_folder + "/" + fov + "_*" + channels[i])
                print(fov_chs[i])
            try:
                assert all(
                    [len(item) for item in fov_chs]
                ), "channels incomplete {}".format(fov)
                print("*" * 40)
                print(fov)
                for i, channel_img in enumerate(fov_chs):
                    fov_chs[i] = proj_channel(channel_img, img_dtype)
                    save_path = os.path.join(plate_folder, fov + channels[i])
                    # print(save_path)
                    imsave(save_path, fov_chs[i])
            except AssertionError:
                continue


if __name__ == "__main__":
    # original datapath
    #acquired_data_path = "/data/trang/covid19_data_CZ8746"
    #formatted_data_folder = "/data/trang/covid19_data_CZ8746_formatted2"
    #data_base = "/data/trang/covid19_data_CZ8746_max_projection2"
    
    # generate csvs, including all-tasks.csv, all-well-tasks.csv
    # where to host the csvs
    #annotation_folder = "/data/trang/covid19_data_CZ8746_annotation2"

    ACQUIRED_DATA_PATH = "/data/trang/211111_COVID19_repurposing_Marianna"
    FORMATTED_DATA_PATH = "/data/trang/211111_COVID19_repurposing_Marianna_formatted"
    IMG_FOLDER = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection"
    ANNOTATION_FOLDER = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection_annotation"
    channels = ["_blue.tiff", "_red.tiff", "_green.tiff", "_yellow.tiff"]
    format_covid_data(ACQUIRED_DATA_PATH, FORMATTED_DATA_PATH)
    generate_max_projection_images(FORMATTED_DATA_PATH, IMG_FOLDER)
    if not os.path.isdir(ANNOTATION_FOLDER):
        os.makedirs(ANNOTATION_FOLDER)
    generate_task_lists(ANNOTATION_FOLDER, IMG_FOLDER)
