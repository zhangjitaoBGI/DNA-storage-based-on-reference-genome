# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 23:36:55 2020

@author: Jitao Zhang
"""
import cv2
import numpy as np
import math
import os
from openpyxl import load_workbook,Workbook

def psnr(img1, img2):
   mse = np.mean( (img1/255. - img2/255.) ** 2 )
   if mse < 1.0e-10:
      return 100
   PIXEL_MAX = 1
   return 20 * math.log10(PIXEL_MAX / math.sqrt(mse))

def calculate_psnr(img1, img2):
    # img1 and img2 have range [0, 255]
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    mse = np.mean((img1 - img2)**2)
    if mse == 0:
        return float('inf')
    return 20 * math.log10(255.0 / math.sqrt(mse))

def ssim(img1, img2):
    C1 = (0.01 * 255)**2
    C2 = (0.03 * 255)**2

    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    kernel = cv2.getGaussianKernel(11, 1.5)
    window = np.outer(kernel, kernel.transpose())

    mu1 = cv2.filter2D(img1, -1, window)[5:-5, 5:-5]  # valid
    mu2 = cv2.filter2D(img2, -1, window)[5:-5, 5:-5]
    mu1_sq = mu1**2
    mu2_sq = mu2**2
    mu1_mu2 = mu1 * mu2
    sigma1_sq = cv2.filter2D(img1**2, -1, window)[5:-5, 5:-5] - mu1_sq
    sigma2_sq = cv2.filter2D(img2**2, -1, window)[5:-5, 5:-5] - mu2_sq
    sigma12 = cv2.filter2D(img1 * img2, -1, window)[5:-5, 5:-5] - mu1_mu2

    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / ((mu1_sq + mu2_sq + C1) *
                                                            (sigma1_sq + sigma2_sq + C2))
    return ssim_map.mean()

def calculate_ssim(img1, img2):
   '''calculate SSIM
   the same outputs as MATLAB's
   img1, img2: [0, 255]
   '''
   if not img1.shape == img2.shape:
       raise ValueError('Input images must have the same dimensions.')
   if img1.ndim == 2:
       return ssim(img1, img2)
   elif img1.ndim == 3:
       if img1.shape[2] == 3:
           ssims = []
           for i in range(3):
               ssims.append(ssim(img1, img2))
           return np.array(ssims).mean()
       elif img1.shape[2] == 1:
           return ssim(np.squeeze(img1), np.squeeze(img2))
   else:
       raise ValueError('Wrong input image dimensions.')
  
'''
img_lena_bw = cv2.imread(r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\small-256\lena-BW-256.bmp")

img_lena_color_256 = cv2.imread(r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\small-256\lena-color-24-256.bmp")

img_set14_lena_512 = cv2.imread(r"E:\work in CNGB\DNA_storage\mmsr\datasets\val_set14\Set14\lenna.png")
'''

'''
img_lena_bw_4 = cv2.imread("E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\out_lena256_bw_4.bmp")
img_lena_bw_2 = cv2.imread("E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\out_lena256_bw_2.bmp")

img_lena_color_222 = cv2.imread("E:/work in CNGB/interim report/image/256/lena-color-256-222.bmp")
img_lena_color_244 = cv2.imread(r"E:\work in CNGB\interim report\image\256\lena-color-256-244.bmp")
img_lena_color_444 = cv2.imread(r"E:\work in CNGB\interim report\image\256\lena-color-256-444.bmp")
'''
'''
img_lena_color_size_6 = cv2.imread(r"E:\work in CNGB\interim report\image\512\size-test\lena-512-244-size_6-seed_2020.bmp")
img_lena_color_size_7 = cv2.imread(r"../size-test/image-result/512/lena-512-244-size_7-seed_2020.bmp")
img_lena_color_size_8 = cv2.imread(r"../size-test/image-result/512/lena-512-244-size_8-seed_2020.bmp")
img_lena_color_size_9 = cv2.imread(r"../size-test/image-result/512/lena-512-244-size_9-seed_2020.bmp")
img_lena_color_size_10 = cv2.imread(r"../size-test/image-result/512/lena-512-244-size_10-seed_2020.bmp")
img_lena_color_size_11 = cv2.imread("../size-test/image-result/512/lena-512-244-size_11-seed_2020.bmp")
img_lena_color_genome = cv2.imread("../size-test/image-result/512/lena-512-244.bmp")
'''
'''
denoisze:
img_lena_color_244_out = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\lenna-244.png")
img_lena_color_244_pred = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\lenna-244-pred.png")
img_monarch = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\monarch.png")
img_monarch_out = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\monarch-244.png")
img_monarch_pred = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\monarch-244-pred.png")
img_ppt3 = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\ppt3.png")
img_ppt3_out = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\ppt3-244.png")
img_ppt3_pred = cv2.imread(r"E:\work in CNGB\interim report\image\512\denoise\ppt3-244-pred.png")
'''

if __name__ == '__main__':    
    original_image_list = os.listdir("./original")
    for image in original_image_list:
        img1 = cv2.imread("./original"+"/"+image)
        img2 = cv2.imread("./genome_based_simulation" + "/"+ "{}_0001.png".format(image[:-4]))
        img3 = cv2.imread("./jpeg_simulation" + "/"+ "{}_0001.jpg".format(image[:-4]))            
    
        print("The result of genome_based image {} is:".format(image[:-4]))
        print(calculate_psnr(img1,img2),calculate_ssim(img1,img2))
        print("The result of jpeg image {} is:".format(image[:-4]))
        print(calculate_psnr(img1,img3),calculate_ssim(img1,img3))
  


