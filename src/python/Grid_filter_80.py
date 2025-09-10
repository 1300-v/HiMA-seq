#import module
import argparse

import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt
import os
import matplotlib
matplotlib.use('Agg')

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--sampleId", required=True, help="sampleId")
args = vars(ap.parse_args())

sampleId = args["sampleId"]

img = cv.imread(sampleId +'.jpg')
image = cv.resize(img,(1500,1500))
# cv_show(image)
image_gray = cv.cvtColor(image,cv.COLOR_BGR2GRAY)

# Gaussian filtering smooths the image
blur = cv.GaussianBlur(image_gray,(5,5),0)
if ((image == 255).sum() >= image.size / 2):
    thresh_img = cv.threshold(image_gray,210,255,cv.THRESH_BINARY)[1]
else:
    thresh_img = cv.threshold(image_gray,200,255,cv.THRESH_BINARY_INV)[1]

# Zoom lossless to 80x80
image_80 = cv.resize(thresh_img,(80,80),interpolation=cv.INTER_LINEAR)
white_pixels = np.where(image_80 == 255)

data = open("position_{}.txt".format(sampleId),'w',encoding="utf-8")
useful_pixel =""
for y, x in zip(white_pixels[0], white_pixels[1]):
    useful_pixel=useful_pixel+","+str(x+1)+"x"+str(y+1)
print(useful_pixel,file=data)
data.close()

plt.figure(figsize=(10,8))
plt.subplot(1,3,1)
plt.imshow(image_gray, "gray")
plt.title("img_gray")
plt.subplot(1,3,2)
plt.imshow(thresh_img,"gray")
plt.title("thresh_img")
plt.subplot(1,3,3)
plt.imshow(image_80, "gray")
plt.title("bin")

plt.savefig(sampleId +'_filtered.jpg')
plt.show()


