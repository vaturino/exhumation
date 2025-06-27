#!/usr/bin/env python3

import cv2
import numpy as np
import argparse
import os
import sys
from tqdm import tqdm
import re

# Function to extract numeric parts for natural sorting
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 script.py <model_name> <video_name>")
        sys.exit(1)

    mod_name = str(sys.argv[1])
    vid_name = str(sys.argv[2])

    image_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/{vid_name}'
    out_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/videos/'
    # image_folder = f"/home/vturino/PhD/projects/exhumation/plots/comparison/{mod_name}"
    # out_folder = f"/home/vturino/PhD/projects/exhumation/plots/comparison/videos/{vid_name}"

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # List and sort PNG images
    image_files = sorted(
        [f for f in os.listdir(image_folder) if f.endswith('.png')],
        key=natural_sort_key
    )

    if not image_files:
        print("No PNG images found in folder.")
        sys.exit(1)

    # Get maximum width/height to standardize video frame size
    max_height, max_width = 0, 0
    for img_file in image_files:
        img = cv2.imread(os.path.join(image_folder, img_file))
        if img is None:
            continue
        h, w, _ = img.shape
        max_height = max(max_height, h)
        max_width = max(max_width, w)

    # Ensure dimensions are even (required for some codecs like H.264)
    if max_width % 2 != 0:
        max_width += 1
    if max_height % 2 != 0:
        max_height += 1

    fps = 3
    output_path = os.path.join(out_folder, f"{vid_name}.mp4")

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # 'avc1' or 'H264' may work too
    video_writer = cv2.VideoWriter(output_path, fourcc, fps, (max_width, max_height))

    for img_file in tqdm(image_files, desc="Writing video"):
        img_path = os.path.join(image_folder, img_file)
        img = cv2.imread(img_path)

        if img is None:
            print(f"Warning: Could not read image {img_path}. Skipping.")
            continue

        h, w, _ = img.shape
        # Pad image to be centered within the frame
        top = (max_height - h) // 2
        left = (max_width - w) // 2
        padded = np.zeros((max_height, max_width, 3), dtype=np.uint8)
        padded[top:top+h, left:left+w] = img


        video_writer.write(padded)

    video_writer.release()
    print(f"✅ Video saved to: {output_path}")

if __name__ == "__main__":
    main()






# #! /usr/bin/python3
# import cv2
# import numpy as np
# import skvideo.io
# import argparse
# import os, sys
# from tqdm import tqdm
# import re

# # Function to extract numeric parts for natural sorting
# def natural_sort_key(s):
#     return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

# def main():
#     mod_name = str(sys.argv[1])
#     vid_name = str(sys.argv[2])

#     image_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/{vid_name}'
#     out_folder = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{mod_name}/videos/'
    
#     if not os.path.exists(out_folder):
#         os.mkdir(out_folder)

#     # List image files and sort them numerically based on the file name
#     image_files = sorted([f for f in os.listdir(image_folder) if f.endswith('.png')], key=natural_sort_key)
#     fnumber = len(image_files)

#     # Determine the maximum dimensions of all images
#     max_height, max_width = 0, 0
#     for img_file in image_files:
#         img = cv2.imread(os.path.join(image_folder, img_file))
#         if img is not None:
#             h, w, _ = img.shape
#             max_height = max(max_height, h)
#             max_width = max(max_width, w)

#     # Create an empty NumPy array with the correct dimensions (max size)
#     out_video = np.empty([fnumber, max_height, max_width, 3], dtype=np.uint8)

#     # Load all images into the NumPy array, padding them to match the largest size
#     for i in tqdm(range(fnumber)):
#         img_path = os.path.join(image_folder, image_files[i])
#         img = cv2.imread(img_path)

#         if img is None:
#             print(f"Error: Image {img_path} could not be loaded.")
#             continue

#         h, w, _ = img.shape
#         # Create a blank image with the max size and fill it with the current image
#         padded_img = np.zeros((max_height, max_width, 3), dtype=np.uint8)
#         padded_img[0:h, 0:w] = img  # Place the original image in the top-left corner
#         out_video[i] = padded_img

#     fps = '2'

#     # Write the output image sequences to a video file
#     skvideo.io.vwrite(f"{out_folder}{vid_name}.mov", out_video, 
#                       inputdict={'-r': fps, "-pix_fmt": "bgr24"},
#                       outputdict={'-f': "mov", "-vcodec": "libx264", "-pix_fmt": "yuv420p", '-vf': "pad=ceil(iw/2)*2:ceil(ih/2)*2"}, 
#                       verbosity=0)

# if __name__ == "__main__":
#     main()
