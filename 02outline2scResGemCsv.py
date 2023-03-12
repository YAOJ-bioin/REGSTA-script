#%%
##Variable definitions
inputPath = r'./data4test/'
imgName = r'img4test.png'
outlineName = r'img4test_cp_outlines'
gemName = r'gem4test'
maxFilterCellArea = "16512.5"
minFilterCellArea = "7.5"
#%%
from cgi import print_form
import os
import sys
import re
import cv2
from PIL import Image
import json
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import colorTools
from colorTools import *
import csv
import pprint
import shutil
from matplotlib import cm
import cv2 as cv
from tqdm import tqdm
#%%
#################################################################################################################
###01 outline2json
#################################################################################################################
def outline2json(file_path,file_name_initial,file_name_final,image_path):
#
    img = Image.open(image_path)
    image_height = str(img.height) + ','
    image_width = img.width
#
    file1 = open(file_path + '/' + file_name_initial, 'r', encoding='utf-8') 
    file2 = open(file_path + '/1.txt', 'w', encoding='utf-8') 
    try:
        for line in file1.readlines():
            if line == '\n':
                line = line.strip("\n")
            file2.write(line)
    finally:
        file1.close()
        file2.close()
#
    ff = open(file_path + '/2.txt','w') 
    with open(file_path + '/1.txt','r') as f:  
        line = f.readlines()  
        for line_list in line:
            line_new ='bioinplant'+line_list
            line_new =line_new.replace('\n','')  
            line_new=line_new+'tnalpnioib'+'\n'  
            str1=r','
            str2=r'\n'
            line_new = re.sub(str1,str2,line_new)
            ff.write(line_new) 
    f.close()
    ff.close()
#
    def joinlns(lns, spliter=","):
        return spliter.join([ln.strip() for ln in lns]) 
    lines = open(file_path + '/2.txt').readlines()
    mergedlines = [joinlns(x, spliter=",") for x in zip(lines[::2], lines[1::2])]
    with open(file_path + '/3.txt', 'w') as handle:
        handle.write("\n".join(mergedlines))
    handle.close()
#
    ff = open(file_path + '/4.txt','w') 
    with open(file_path + '/3.txt','r') as f:  
        line = f.readlines()  
        for line_list in line:
            line_new ='['+line_list
            line_new =line_new.replace('\n','')  
            line_new=line_new+'],'+'\n'  
            str1=r'\[bioinplant'
            str2=r'{\n"label": "cell",\n"points": [\n['
            line_new=re.sub(str1,str2,line_new)
            str3=r'tnalpnioib],'
            str4=r']\n],\n"group_id": null,\n"shape_type": "polygon",\n"flags": {}\n},'
            line_new=re.sub(str3,str4,line_new)
            ff.write(line_new) 
    f.close()
    ff.close()
#
    ff = open(file_path + '/5.txt','w+') 
    with open(file_path + '/4.txt', 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        ff.write('{\n"version": "5.0.1",\n"flags": {},\n"shapes": [\n'+content+'],\n"imagePath":"' + str(image_path) + '",\n"imageHeight":' + str(image_height) + '\n"imageWidth":' + str(image_width) + '\n}')
    ff.close()
#
    ff = open(file_path + '/' + file_name_final,'w') 
    with open(file_path + '/5.txt','r') as f:  
        line = f.read()  
        str1=r'},\n],'
        str2=r'}\n],'
        line_new=re.sub(str1,str2,line)
        ff.write(line_new) 
    f.close()
    ff.close()
#
    for num in range(1,6):
        os.remove(file_path + '/' + str(num) + '.txt')

if __name__ == "__main__":
    file_path = inputPath
    file_name_initial = outlineName + r'.txt'
    file_name_final = outlineName + r'.json'
    image_path = inputPath + imgName
    outline2json(file_path,file_name_initial,file_name_final,image_path)

#################################################################################################################
###02 json2coloredJson
#################################################################################################################
import json
import random
def gen_color_df(color_list):
    pd_str_list = []
    for i, x in enumerate(color_list):      
        pd_str_list.append(str(x))         

    df_color = pd.concat([pd.DataFrame(color_list, columns=['r', 'g', 'b']),     
                          pd.DataFrame(pd_str_list, columns=['str'])], axis=1)   
    df_color.to_csv('color_list_int.csv')  
    return


def get_color_list(num):
    color_list = []

    for i in range(1, num+1):
        rgb = int2rgb(i)
        color_list.append(rgb)
    # 不shuffle 按照顺序染色
    # random.seed(0)
    # random.shuffle(color_list)
    return color_list

def get_color_list_onlyshow(num):
    color_list = []
    color_camp = plt.cm.get_cmap('hsv', num)
    for i in range(num):
        tmp_camp = color_camp(i)
        pretty_camp = (tmp_camp[0] * 255, tmp_camp[1] * 255, tmp_camp[2] * 255)
        color_list.append(pretty_camp)
    random.seed(0)
    random.shuffle(color_list)
    return color_list


def convertPolygonToColoredMask(jsonfilePath, scaling=1):
    with open(jsonfilePath, "r", encoding='utf-8') as jsonf:
        jsonData = json.load(jsonf)
        img_h = jsonData["imageHeight"] * scaling
        img_w = jsonData["imageWidth"] * scaling
        mask = np.zeros((img_h, img_w, 3), np.uint8)

        num_sum = len(jsonData["shapes"])  # 图片中目标的数量****4965******
        color_list = get_color_list(num_sum)
        gen_color_df(color_list)  # 生成csv保存

        i = 0
        for obj in jsonData["shapes"]:
            label = obj["label"]
            polygonPoints = obj["points"]
            polygonPoints = np.array(polygonPoints, np.int32) * scaling
            # print("+" * 50, "\n", polygonPoints)
            # print(label)
            cv2.drawContours(mask, [polygonPoints], -1, color_list[i], -1)
            i += 1
    maskRGB = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
    return maskRGB


def convertPolygonToColoredJson(jsonfilePath):
    with open(jsonfilePath, "r", encoding='utf-8') as jsonf:
        jsonData = json.load(jsonf)
        img_h = jsonData["imageHeight"]
        img_w = jsonData["imageWidth"]
        mask = np.zeros((img_h, img_w, 3), np.uint8)

        num_sum = len(jsonData["shapes"])  # 图片中目标的数量
        color_list = get_color_list(num_sum)
        # gen_color_df(color_list)  # 生成csv保存

        i = 0
        for obj in jsonData["shapes"]:
            label = obj["label"]
            polygonPoints = obj["points"]
            polygonPoints = np.array(polygonPoints, np.int32)

            cv2.drawContours(mask, [polygonPoints], -1, color_list[i], -1)
            obj["preCellRGB"] = color_list[i]
            obj["preCellId"] = colorTools.rgb2int(color_list[i][0], color_list[i][1], color_list[i][2])

            # 计算面积周长也放在json里
            M = cv2.moments(polygonPoints)
            gtCellArea = M['m00']
            gtCellPerimeter = cv2.arcLength(polygonPoints, True)

            obj['preCellArea'] = gtCellArea
            obj['preCellPerimeter'] = gtCellPerimeter

            i += 1

    maskRGB = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
    return maskRGB, jsonData


if __name__ == "__main__":
    jsonfilePath = inputPath + outlineName + r'.json'
    maskSavePath = inputPath + outlineName + "_colored.png"
    jsonSavePath = inputPath + outlineName + "_colored.json"

    mask, jsonData = convertPolygonToColoredJson(jsonfilePath)

    cv2.imwrite(maskSavePath, mask)
    cv2.destroyAllWindows()

    with open(jsonSavePath, 'w', encoding='utf-8') as fp:
        json.dump(jsonData, fp)
    os.remove(inputPath + outlineName + '.json')

#################################################################################################################
###03 coloredJson2coloredNewJson
#################################################################################################################
def json2csv(input,output_json,output_csv):

    os.chdir(inputPath)
#
    f1 = open(input,'r+')
    f2 = open('1.json','w+')
    str1=r'{"label"'
    str2=r'\n{"label"'
    for ss in f1.readlines():
        tt=re.sub(str1,str2,ss)
        f2.write(tt)
    f1.close()
    f2.close()
#
    f1 = open('1.json','r+')
    f2 = open('2.json','w+')
    str1=r', "imagePath":'
    str2=r'\n, "imagePath":'
    for ss in f1.readlines():
        tt=re.sub(str1,str2,ss)
        f2.write(tt)
    f1.close()
    f2.close()
#
    f1 = open('2.json','r+')
    f2 = open('3.json','w+')
    str1=r'{"version": "5.0.1", "flags": {}, "shapes": '
    str2=r''
    for ss in f1.readlines():
        tt=re.sub(str1,str2,ss)
        f2.write(tt)
    f1.close()
    f2.close()
#
    with open('3.json','r') as r:
        lines=r.readlines()
    with open(output_json,'w') as w:
        for l in lines:
            if 'imagePath' not in l:
                w.write(l) 
#
    for num in range(1,4):
        os.remove(str(num) + '.json')
#######
#
    with open(output_json,'r',encoding='utf8')as fp:
        json_data = json.load(fp)
    csv_file = open(output_csv, 'w')
    sheet_title = json_data[0].keys()
    json_values = []
    for dict in json_data:
        json_values.append(dict.values())
    csv_writer = csv.writer(csv_file)

    # 4.2 写入表头
    csv_writer.writerow(sheet_title)
    # 4.3 写入内容
    csv_writer.writerows(json_values)

    # 5.关闭文件
    csv_file.close()
#    json_file.close()


    data = pd.read_csv(output_csv)
    res = data.dropna(how="all")
    res.to_csv(output_csv, index=False)
    #os.remove(output_json)


if __name__ == "__main__":
    
    input = inputPath + outlineName + "_colored.json"
    output_json = inputPath + outlineName + "_colored_new.json"
    output_csv = inputPath + outlineName + "_colored.csv"
    json2csv(input,output_json,output_csv)
    os.remove(inputPath + outlineName + '_colored.png')
    os.remove(inputPath + outlineName + '_colored.json')
#################################################################################################################
###04 coloredNewJson2newTxt
#################################################################################################################
os.chdir(inputPath)

with open(inputPath + outlineName + "_colored_new.json",'r',encoding='utf8')as fp:
    json_data = json.load(fp)
json=pd.DataFrame(json_data)
#
df = json[(json['preCellArea'] <= float(maxFilterCellArea) )&(json['preCellArea'] >= float(minFilterCellArea) )]
dtf = df['points']
dtf.to_csv('1.txt', sep='\t', index=False)
#
ff = open("2.txt",'w') 
with open('1.txt','r') as f:  
    line = f.read()  
    str1=r'\['
    str2=r''
    line_new=re.sub(str1,str2,line)
    ff.write(line_new) 
f.close()
ff.close()
#
ff = open("3.txt",'w') 
with open('2.txt','r') as f:  
    line = f.read()  
    str1=r'\]'
    str2=r''
    line_new=re.sub(str1,str2,line)
    ff.write(line_new) 
f.close()
ff.close()
#
ff = open("4.txt",'w') 
with open('3.txt','r') as f:  
    line = f.read()  
    str1=r' '
    str2=r''
    line_new=re.sub(str1,str2,line)
    ff.write(line_new) 
f.close()
ff.close()
#
ff = open(inputPath + outlineName + "_colored_new.txt",'w') 
with open('4.txt','r') as f:  
    line = f.read()  
    str1=r'points\n'
    str2=r''
    line_new=re.sub(str1,str2,line)
    ff.write(line_new) 
f.close()
ff.close()
#
for num in range(1,5):
    os.remove(inputPath + str(num) + '.txt')
os.remove(inputPath + outlineName + '_colored_new.json')
os.remove(inputPath + outlineName + '_colored.csv')

#################################################################################################################
###05 newTxt2finalJson
#################################################################################################################
if __name__ == "__main__":
    file_path = inputPath
    file_name_initial = outlineName + r'_colored_new.txt'
    file_name_final = outlineName + r'4align.json'
    image_path = inputPath + imgName
    outline2json(file_path,file_name_initial,file_name_final,image_path)
    os.remove(inputPath + outlineName + '_colored_new.txt')
#    source = outlineName + r'4align.json'
#    target = savePath
#    shutil.copy(source, target)

#################################################################################################################
###06.1 finalJson2handrollMask
#################################################################################################################
import json
import random
def convertPolygonToMask(jsonfilePath, scaling=1):
    with open(jsonfilePath, "r", encoding='utf-8') as jsonf:
        jsonData = json.load(jsonf)
        img_h = jsonData["imageHeight"]*scaling
        img_w = jsonData["imageWidth"]*scaling
        mask = np.zeros((img_h, img_w), np.uint8)
        num = 0
        for obj in jsonData["shapes"]:
            label = obj["label"]
            polygonPoints = obj["points"]
            polygonPoints = np.array(polygonPoints, np.int32)*scaling
            num += 1
            cv.drawContours(mask, [polygonPoints], -1, (255), -1)

    return mask


def main():
    jsonfileFolder = r"K:\imageData\colorR\dataset\label"
    maskSaveFolder = r"K:\imageData\colorR\dataset\mask"

    for jsonfile in os.listdir(jsonfileFolder):
        jsonfilePath = os.path.join(jsonfileFolder, jsonfile)
        mask = convertPolygonToMask(jsonfilePath, 4)
        maskName = jsonfile.split(".")[0] + ".png"
        maskPath = os.path.join(maskSaveFolder, maskName)
        cv.imwrite(maskPath, mask)

if __name__ == "__main__":
    # main()
    jsonfilePath = inputPath + outlineName + r'4align.json'
    maskSavePath = inputPath + outlineName + r"_mask_handroll.png"
    mask = convertPolygonToMask(jsonfilePath, 1)

    _, th = cv.threshold(mask, 0, 255, cv.THRESH_BINARY)
    cv.imwrite(maskSavePath, mask)
    cv.destroyAllWindows()

#################################################################################################################
###06.2 finalJson2coloredMaskInt
#################################################################################################################
def gen_color_df(color_list):
    pd_str_list = []
    for i, x in enumerate(color_list):
        pd_str_list.append(str(x))

    df_color = pd.concat([pd.DataFrame(color_list, columns=['r', 'g', 'b']),
                          pd.DataFrame(pd_str_list, columns=['str'])], axis=1)
    df_color.to_csv('color_list_int.csv')
    return


def get_color_list(num):
    color_list = []

    for i in range(num):
        rgb = int2rgb(i)
        color_list.append(rgb)
    return color_list

def get_color_list_onlyshow(num):
    color_list = []
    color_camp = plt.cm.get_cmap('hsv', num)
    for i in range(num):
        tmp_camp = color_camp(i)
        pretty_camp = (tmp_camp[0] * 255, tmp_camp[1] * 255, tmp_camp[2] * 255)
        color_list.append(pretty_camp)
    random.seed(0)
    random.shuffle(color_list)
    return color_list


def convertPolygonToColoredMask(jsonfilePath, scaling=1):
    with open(jsonfilePath, "r", encoding='utf-8') as jsonf:
        jsonData = json.load(jsonf)
        img_h = jsonData["imageHeight"] * scaling
        img_w = jsonData["imageWidth"] * scaling
        mask = np.zeros((img_h, img_w, 3), np.uint8)

        num_sum = len(jsonData["shapes"]) 
        i = 0
        color_list = get_color_list(num_sum)
        gen_color_df(color_list) 
        for obj in jsonData["shapes"]:
            label = obj["label"]
            polygonPoints = obj["points"]
            polygonPoints = np.array(polygonPoints, np.int32) * scaling
            # print("+" * 50, "\n", polygonPoints)
            # print(label)
            cv2.drawContours(mask, [polygonPoints], -1, color_list[i], -1)
            i += 1
    maskRGB = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
    return maskRGB


if __name__ == "__main__":
    # main()
    jsonfilePath = inputPath + outlineName + r'4align.json'
    maskSavePath = inputPath + outlineName + r"_colored_mask_int.png"
    mask = convertPolygonToColoredMask(jsonfilePath, 1) 

    cv2.imwrite(maskSavePath, mask)
    cv2.destroyAllWindows()
    
    os.remove(inputPath + outlineName + r'4align.json')

#################################################################################################################
###07 gem2csv
#################################################################################################################
if __name__ == '__main__':
    with open(inputPath + outlineName + r'.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        for line in open(inputPath + gemName + r'.gem', 'r'):
            line = line
            new_line = line.strip().split()
            f_csv.writerow(new_line)

#################################################################################################################
###08 csv2stImgLocCsv&st4alignImg
#################################################################################################################
if __name__ == '__main__':
    df_st = pd.read_csv(inputPath + outlineName + r'.csv')
    np_xy = df_st.iloc[:, 1:3].values
    np_x = np_xy[:, 0]
    np_y = np_xy[:, 1]
    np_x -= min(np_x)
    np_y -= min(np_y)

    df_xy = pd.concat([pd.DataFrame(np_x, columns=['img_x']), pd.DataFrame(np_y, columns=['img_y'])], axis=1)
    df_new = pd.concat([df_st, df_xy], axis=1)
    df_new.to_csv(inputPath + outlineName + '_st_imgloc.csv', index=False)

    np_mat = np.zeros((max(np_y) + 1, max(np_x) + 1))
    for i in range(len(np_x)):
        np_mat[np_y[i], np_x[i]] = 1
        #print(np_y[i], np_x[i])
        #print()
    np_mat *= 255

    image = Image.fromarray(np_mat)
    image = image.convert('L')
    image.save(inputPath + outlineName + '_st4align.png')

#################################################################################################################
###09 ____2mask_resize2st_colored
#################################################################################################################

def order_points(pts):
    ''' sort rectangle points by clockwise '''
    sort_x = pts[np.argsort(pts[:, 0]), :]
    Left = sort_x[:2, :]
    Right = sort_x[2:, :]
    Left = Left[np.argsort(Left[:, 1])[::-1], :]
    Right = Right[np.argsort(Right[:, 1]), :]
    return np.concatenate((Left, Right), axis=0)


def get_rectangle_pts(res):
    contours, hierarchy = cv2.findContours(res, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    area_lst = []
    for i, cnt in enumerate(contours):
        area = cv2.contourArea(cnt)
        area_lst.append(area)
    maxi = np.argmax(np.array(area_lst))
    maxContour = contours[maxi]
    x, y, w, h = cv2.boundingRect(maxContour)
    return x, y, w, h, maxContour


def show_img_now(img, title):
    plt.imshow(img)
    plt.title(title)
    plt.show()

if __name__ == '__main__':
    mask_img = cv2.imread(inputPath + outlineName + r"_mask_handroll.png")
    mask_img_gray = cv2.cvtColor(mask_img, cv2.COLOR_BGR2GRAY)

    mask_img = cv2.imread(inputPath + outlineName + r"_colored_mask_int.png")
    mask_img_rgb = cv2.cvtColor(mask_img, cv2.COLOR_BGR2RGB)

    _, mask_res = cv2.threshold(mask_img_gray, 150, 255, cv2.THRESH_BINARY)
    kernel = np.ones((55, 55), np.uint8)  
    mask_res = cv2.morphologyEx(mask_res, cv2.MORPH_CLOSE, kernel)

    x, y, w, h, maxContour = get_rectangle_pts(mask_res)
    mask_with_rectangle = cv2.rectangle(mask_img, (x, y), (x + w, y + h), (255, 255, 255), 20)
    print('mask img real rectangle shape:=', x, y, w, h, '(xywh)')
    colored_crop = mask_img_rgb[y:y + h, x:x + w]

    st_img = cv2.imread(inputPath + outlineName + '_st4align.png')
    mask_resize2st = cv2.resize(colored_crop, (st_img.shape[1], st_img.shape[0]))

    syn_img = cv2.addWeighted(mask_resize2st, 0.2, st_img, 0.8, 0)

    syn_img = cv2.cvtColor(syn_img, cv2.COLOR_BGR2RGB)
    mask_resize2st = cv2.cvtColor(mask_resize2st, cv2.COLOR_BGR2RGB)
    cv2.imwrite(inputPath + outlineName + '_syn_img_colored.png', syn_img)
    cv2.imwrite(inputPath + outlineName + '_mask_resize2st_colored.png', mask_resize2st)

    os.remove(inputPath + outlineName + r"_mask_handroll.png")
    os.remove(inputPath + outlineName + r"_colored_mask_int.png")
    os.remove(inputPath + outlineName + r'.csv')


#################################################################################################################
###10 ___2scResGemCsv
#################################################################################################################

if __name__ == '__main__':
    colored_mask = cv2.imread(inputPath + outlineName + '_mask_resize2st_colored.png')
    colored_mask = cv2.cvtColor(colored_mask, cv2.COLOR_BGR2RGB)

    st_img = cv2.imread(inputPath + outlineName + '_st4align.png')

    img_size = st_img.shape
    img_h0 = img_size[0]  # 3679
    img_w0 = img_size[1]  # 2470

    df_st = pd.read_csv(inputPath + outlineName + '_st_imgloc.csv')

    cell_id_list = []
    for i in tqdm(range(df_st.shape[0])):
        df_row_tmp = df_st.iloc[i]
        y_tmp = df_row_tmp.loc['img_y']
        x_tmp = df_row_tmp.loc['img_x']

        mask_value = colored_mask[y_tmp][x_tmp]
        cell_id_tmp = colorTools.rgb2int(mask_value[0], mask_value[1], mask_value[2])

        cell_id_list.append(cell_id_tmp)

    df_cell_id = pd.DataFrame(cell_id_list, columns=['cell_id'])
    cell_id_res = pd.concat([df_st, df_cell_id], axis=1)

    cell_id_res.to_csv(inputPath + outlineName +'_scResGem.csv')

    print()

    os.remove(inputPath + outlineName + '_mask_resize2st_colored.png')
    os.remove(inputPath + outlineName + '_st4align.png')
    os.remove(inputPath + outlineName + '_st_imgloc.csv')
    os.remove(inputPath + outlineName + '_syn_img_colored.png')
    os.remove(inputPath + 'color_list_int.csv')