# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# ### This is to format the acquired image plates to the format of covid-annotation and then generate the all-tasks.csv and all-well-tasks.csv
# ##### the csv files should be copied to the covid-annotation folder

import os
from shutil import copyfile


def format_covid_data(data_path, final_data_folder):
    os.makedirs(final_data_folder, exist_ok=True)
    plates = list(
        filter(
            lambda item: item if not item.startswith(".") else None,
            os.listdir(data_path),
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
        "ch1": "blue.tiff",
        "ch2": "red.tiff",
        "ch3": "green.tiff",
        "ch4": "yellow.tiff",
    }
    for plate in plates:
        plate_folder = os.path.join(final_data_folder, plate)
        os.makedirs(plate_folder, exist_ok=True)
        # plate_images[plate] = all_images
        all_images = []
        plate_path = os.path.join(data_path, plate)
        image_folder = os.path.join(plate_path, "Images")

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
        fov_index = 1
        for i, item in enumerate(fovs):
            row_index = row_mapping[item[:3]]
            column_index = str(int(item[4:6]))
            if not i:
                old_row = row_index
                old_column = column_index
            else:
                if not (
                    (old_row == row_index) and (old_column == column_index)
                ):
                    # print('Nej')
                    old_row = row_index
                    old_column = column_index
                    fov_index = 1
            mapping_dict[item] = "".join(
                [
                    row_index,
                    column_index,
                    "_",
                    str(fov_index).zfill(FOV_INDEX_LENGTH),
                ]
            )
            # print(''.join([row_index, column_index, '_', str(fov_index)]))
            fov_index += 1
        for image in all_images:
            source_path = os.path.join(data_path, plate, "Images", image)
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


if __name__ == "__main__":
    # original datapath
    data_path = "/data/trang/covid19_data/"
    data_base = "/data/trang/formatted_covid_images"
    # generate csvs, including all-tasks.csv, all-well-tasks.csv
    # where to host the csvs
    csv_folder = "/data/trang/test_removable"
    # set the fov index length
    # not fov_index_length should be set depending on how many fovs in one
    # well. like 9 fov and 3 z-slice image, 3*9=27, then fov_index_length = 2.
    # this is to ensure the z-stack be together in the annotation tool
    FOV_INDEX_LENGTH = 2
    format_covid_data(data_path, data_base)
    generate_task_lists(csv_folder, data_base)
