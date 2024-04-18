import os
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import lindhardt as lndh   
import cv2
import pytesseract

plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"xtick.labelsize": 15})
plt.rcParams.update({"ytick.labelsize": 15})
plt.rcParams.update({"axes.labelsize": 15})
def extract_values_from_gui(root_folder):
    Energy = []
    Recoil_phonon = []
    Ion_phonon = []
    Recoil_vacancy = []
    Ion_vacancy = []
    Recoil_ionization = []
    Ion_ionization = []
    Vacancy = []
    for item in sorted(os.listdir(root_folder), key=lambda x: (float(re.search(r'\d+\.\d+|\d+', x).group()) if re.search(r'\d+\.\d+|\d+', x) else 0, x)):
        item_path = os.path.join(root_folder, item)
        if item == ".vscode":
            continue
        elif item == "GUI_summary.txt":
            continue
        energy =float(item.split('_')[0])
        print(energy)
        Energy.append(energy)
        if os.path.isdir(item_path):
            # Iterate over all the files in the subfolder
            for file_name in os.listdir(item_path):         
                if (file_name==str(energy)+".png" or file_name ==str(energy).replace(".", "_")+".png"):
                    print(file_name)
                    img = cv2.imread(fr'{item_path}/{file_name}')
                    img_crop_vacancies = img[165:190, 1300:1385]
                    img_crop = img[370:450, 1300:1385]
                    for idx, img in enumerate([img_crop_vacancies, img_crop]):
                        img_hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
                        lower_red = np.array([0,50,50])
                        upper_red = np.array([10,255,255])
                        mask0 = cv2.inRange(img_hsv, lower_red, upper_red)

                        # upper mask (170-180)
                        lower_red = np.array([170,50,50])
                        upper_red = np.array([180,255,255])
                        mask1 = cv2.inRange(img_hsv, lower_red, upper_red)

                        # join my masks
                        mask = mask0+mask1
                        # Use pytesseract to extract the numbers from the image
                        kernel = np.ones((2,2), np.uint8)
                        #mask = cv2.dilate(mask, kernel, iterations=1)
                        rect_kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
                        dilation = cv2.dilate(mask, rect_kernel, iterations = 3)
                        # cv2.imshow("mask",   mask)
                        # cv2.imshow('dilation', dilation)
                        # cv2.waitKey(0)
                        # exit()
                        contours, hierarchy = cv2.findContours(dilation, cv2.RETR_EXTERNAL,
                                                cv2.CHAIN_APPROX_NONE)

                        im2 = mask.copy()
                        i=0
                        for cnt in contours:
                            x, y, w, h = cv2.boundingRect(cnt)                        
                            # Draw the bounding box on the text area
                            rect=cv2.rectangle(im2, (x, y), (x + w, y + h), (0, 255, 0), 2)
                            # Crop the bounding box area
                            scale_factor = 4
                            cropped = mask[y:y + h, x:x + w]
                            cropped = cv2.bitwise_not(cropped)
                            cv2.imwrite("picture"+str(i)+".png",cropped)  
                            image = cv2.imread("picture"+str(i)+".png")
                            image = cv2.resize(image, (0,0), fx=4, fy=4)
                            #image = imutils.resize(image, width=300)
                            thresh = cv2.threshold(image, 150, 255, cv2.THRESH_BINARY_INV)[1]

                            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3))
                            close = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)

                            result = 255 - close 
                            result = cv2.GaussianBlur(result, (5,5), 0)

                        
                            #cv2.imwrite('rectanglebox.jpg',rect)
                            config='--psm 13 digits' 
                            #config='--psm 6 --oem 3 -c tessedit_char_whitelist=0123456789.' 
                            text = pytesseract.image_to_string(result, config = config)            
                            # Adding the text to the file
                            if idx==1:
                                if i==0:
                                    recoil_phonon = float(text)
                                    Recoil_phonon.append(recoil_phonon)
                                elif i==1:
                                    ion_phonon = float(text)
                                    Ion_phonon.append(ion_phonon)
                                elif i==2:
                                    recoil_vacancy = float(text)
                                    Recoil_vacancy.append(recoil_vacancy)
                                elif i==3:
                                    ion_vacancy = float(text)
                                    Ion_vacancy.append(ion_vacancy)
                                elif i==4:
                                    recoil_ionization = float(text)
                                    Recoil_ionization.append(recoil_ionization)
                                elif i==5:
                                    ion_ionization = float(text)
                                    Ion_ionization.append(ion_ionization)
                            elif idx==0:
                                if i==0:
                                    vacancy = float(text)
                                    Vacancy.append(vacancy)
                            i+=1
    np.savetxt("GUI_summary.txt", np.column_stack((Energy, Ion_ionization, Recoil_ionization,Ion_vacancy, Recoil_vacancy,  Ion_phonon, Recoil_phonon,  Vacancy)), fmt="%s", delimiter="\t")
                                                 
    return
extract_values_from_gui("Silicon_lattice_4")

