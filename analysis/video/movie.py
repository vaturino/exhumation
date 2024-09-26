#! /usr/bin/python3


# import cv2
# import numpy as np
# import skvideo.io
# import skvideo.datasets
# bbb = skvideo.datasets.bigbuckbunny()
# import argparse
# import json as json
# import os, sys
# from tqdm import tqdm



# def main():
#     mod_name = str(sys.argv[1])
#     vid_name = str(sys.argv[2])

#     json_loc='/home/vturino/PhD/projects/exhumation/pyInput/'


#     image_folder=f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/{vid_name}'
#     out_folder=f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/videos/'
#     if not os.path.exists(out_folder):
#             os.mkdir(out_folder)
    
#     fnumber = len(os.listdir(image_folder))


#     out_video =  np.empty([fnumber, 9496,12354,3], dtype = np.uint8)
#     out_video =  out_video.astype(np.uint8)

#     for i in tqdm(range(1, fnumber)):
#         img = cv2.imread(f"{image_folder}/{i}" + '.png')
#         # print(f"{image_folder}/{i}" + '.png')
#         out_video[i] = img

#     fps = str('2')
#     # Writes the the output image sequences in a video file
#     # skvideo.io.vwrite(f"{out_folder}{vid_name}.mov", out_video, inputdict={'-r': fps, "-pix_fmt": "bgr24"}, outputdict={'-f': "mov", "-vcodec": "libx264", "-pix_fmt": "yuv420p", '-vf': "pad=ceil(iw/2)*2:ceil(ih/2)*2"}, verbosity=0)
#     skvideo.io.vwrite(f"{out_folder}{vid_name}.mov", out_video, inputdict={'-r': fps, "-pix_fmt": "bgr24"}, outputdict={'-f': "mov", "-vcodec": "libx264", "-pix_fmt": "yuv420p"}, verbosity=0)


#     # # Closes the video file
#     # skvideo.io.vclose()

# if __name__ == "__main__":
#     main()


#! /usr/bin/python3
import cv2
import numpy as np
import skvideo.io
import argparse
import os, sys
from tqdm import tqdm
import re

# Function to extract numeric parts for natural sorting
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

def main():
    mod_name = str(sys.argv[1])
    vid_name = str(sys.argv[2])

    image_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/{vid_name}'
    out_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/videos/'
    
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    # List image files and sort them numerically based on the file name
    image_files = sorted([f for f in os.listdir(image_folder) if f.endswith('.png')], key=natural_sort_key)
    fnumber = len(image_files)

    # Determine the maximum dimensions of all images
    max_height, max_width = 0, 0
    for img_file in image_files:
        img = cv2.imread(os.path.join(image_folder, img_file))
        if img is not None:
            h, w, _ = img.shape
            max_height = max(max_height, h)
            max_width = max(max_width, w)

    # Create an empty NumPy array with the correct dimensions (max size)
    out_video = np.empty([fnumber, max_height, max_width, 3], dtype=np.uint8)

    # Load all images into the NumPy array, padding them to match the largest size
    for i in tqdm(range(fnumber)):
        img_path = os.path.join(image_folder, image_files[i])
        img = cv2.imread(img_path)

        if img is None:
            print(f"Error: Image {img_path} could not be loaded.")
            continue

        h, w, _ = img.shape
        # Create a blank image with the max size and fill it with the current image
        padded_img = np.zeros((max_height, max_width, 3), dtype=np.uint8)
        padded_img[0:h, 0:w] = img  # Place the original image in the top-left corner
        out_video[i] = padded_img

    fps = '2'

    # Write the output image sequences to a video file
    skvideo.io.vwrite(f"{out_folder}{vid_name}.mov", out_video, 
                      inputdict={'-r': fps, "-pix_fmt": "bgr24"},
                      outputdict={'-f': "mov", "-vcodec": "libx264", "-pix_fmt": "yuv420p", '-vf': "pad=ceil(iw/2)*2:ceil(ih/2)*2"}, 
                      verbosity=0)

if __name__ == "__main__":
    main()
